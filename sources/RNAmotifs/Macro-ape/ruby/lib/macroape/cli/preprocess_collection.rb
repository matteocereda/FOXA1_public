require_relative '../../macroape'
require 'yaml'
require 'shellwords'

module Macroape
  module CLI
    module PreprocessCollection

      def self.main(argv)
        doc = <<-EOS.strip_doc
          Command-line format:
          #{run_tool_cmd} <file or folder with PWMs or .stdin with filenames> <output file> [options]

          Options:
            [-p <list of P-values>] - comma separated(no spaces allowed) list of P-values to precalculate thresholds
            [-d <rough discretization>,<precise discretization>] - set discretization rates, comma delimited (no spaces allowed), order doesn't matter
            [--silent] - hide current progress information during scan (printed to stderr by default)
            [--pcm] - treat the input file as Position Count Matrix. PCM-to-PWM transformation to be done internally.
            [--boundary lower|upper] Upper boundary (default) means that the obtained P-value is greater than or equal to the requested P-value
            [-b <background probabilities] ACGT - 4 numbers, comma-delimited(spaces not allowed), sum should be equal to 1, like 0.25,0.24,0.26,0.25

          The tool preprocesses and stores Macroape motif collection in the specified YAML-file.

          Example:
            #{run_tool_cmd} ./motifs  collection.yaml -p 0.001,0.0005,0.0001 -d 1,10 -b 0.2,0.3,0.3,0.2
        EOS

        if argv.empty? || ['-h', '--h', '-help', '--help'].any?{|help_option| argv.include?(help_option)}
          STDERR.puts doc
          exit
        end

        data_model = argv.delete('--pcm') ? Bioinform::PCM : Bioinform::PWM

        default_pvalues = [0.0005]
        background = [1,1,1,1]
        rough_discretization = 1
        precise_discretization = 10
        max_hash_size = 10000000

        data_source = argv.shift
        output_file = argv.shift

        raise 'No input. You should specify file or folder with pwms' unless data_source
        raise "Error! File or folder #{data_source} doesn't exist" unless Dir.exist?(data_source) || File.exist?(data_source) || data_source == '.stdin'
        raise 'You should specify output file'  unless output_file

        pvalues = []
        silent = false
        pvalue_boundary = :upper

        until argv.empty?
          case argv.shift
            when '-b'
              background = argv.shift.split(',').map(&:to_f)
              raise 'background should be symmetric: p(A)=p(T) and p(G) = p(C)' unless background == background.reverse
            when '-p'
              pvalues = argv.shift.split(',').map(&:to_f)
            when '-d'
              rough_discretization, precise_discretization = argv.shift.split(',').map(&:to_f).sort
            when '--max-hash-size'
              max_hash_size = argv.shift.to_i
            when '--silent'
              silent = true
            when '--boundary'
              pvalue_boundary = argv.shift.to_sym
              raise 'boundary should be either lower or upper'  unless  pvalue_boundary == :lower || pvalue_boundary == :upper
            end
        end
        pvalues = default_pvalues  if pvalues.empty?

        collection = Bioinform::Collection.new(rough_discretization: rough_discretization,
                                precise_discretization: precise_discretization,
                                background: background,
                                pvalues: pvalues)

        data_source = data_source.gsub("\\",'/')
        if File.directory?(data_source)
          motifs = Dir.glob(File.join(data_source,'*')).sort.map do |filename|
            pwm = data_model.new(File.read(filename))
            pwm.name ||= File.basename(filename, File.extname(filename))
            pwm
          end
        elsif File.file?(data_source)
          input = File.read(data_source)
          motifs = data_model.split_on_motifs(input)
        elsif data_source == '.stdin'
          filelist = $stdin.read.shellsplit
          motifs = []
          filelist.each do |filename|
            motif = data_model.new(File.read(filename))
            motif.name ||= File.basename(filename, File.extname(filename))
            motifs << motif
          end
        else
          raise "Specified data source `#{data_source}` is neither directory nor file nor even .stdin"
        end

        pwms = motifs.map(&:to_pwm)

        pwms.each_with_index do |pwm,index|
          STDERR.puts "Motif #{pwm.name}, length: #{pwm.length} (#{index+1} of #{pwms.size}, #{index*100/pwms.size}% complete)"  unless silent

          # When support of onefile collections is introduced - then here should be check if name exists.
          # Otherwise it should skip motif and tell you about this
          # Also two command line options to fail on skipping or to skip silently should be included

          info = OpenStruct.new(rough: {}, precise: {})
          pwm.set_parameters(background: background, max_hash_size: max_hash_size)
          skip_motif = false


          fill_rough_infos = ->(pvalue, threshold, real_pvalue) do
            if real_pvalue == 0
              $stderr.puts "#{pwm.name} at pvalue #{pvalue} has threshold that yields real-pvalue 0 in rough mode. Rough calculation will be skipped"
            else
              info.rough[pvalue] = threshold / rough_discretization
            end
          end

          fill_precise_infos = ->(pvalue, threshold, real_pvalue) do
            if real_pvalue == 0
              $stderr.puts "#{pwm.name} at pvalue #{pvalue} has threshold that yields real-pvalue 0 in precise mode. Motif will be excluded from collection"
              skip_motif = true
            else
              info.precise[pvalue] = threshold / precise_discretization
            end
          end

          if pvalue_boundary == :lower
            pwm.discrete(rough_discretization).thresholds(*pvalues, &fill_rough_infos)
          else
            pwm.discrete(rough_discretization).weak_thresholds(*pvalues, &fill_rough_infos)
          end

          if pvalue_boundary == :lower
            pwm.discrete(precise_discretization).thresholds(*pvalues, &fill_precise_infos)
          else
            pwm.discrete(precise_discretization).weak_thresholds(*pvalues,&fill_precise_infos)
          end
          collection.add_pm(pwm, info)  unless skip_motif
        end
        STDERR.puts "100% complete. Saving results"  unless silent
        File.open(output_file, 'w') do |f|
          f.puts(collection.to_yaml)
        end
        puts OutputInformation.new{|infos|
          infos.add_parameter('P', 'P-value list', pvalues.join(','))
          infos.add_parameter('VR', 'discretization value, rough', rough_discretization)
          infos.add_parameter('VP', 'discretization value, precise', precise_discretization)
          infos.add_parameter('PB', 'P-value boundary', pvalue_boundary)
          infos.background_parameter('B', 'background', background)
        }.result
      rescue => err
        STDERR.puts "\n#{err}\n#{err.backtrace.first(5).join("\n")}\n\nUse --help option for help\n\n#{doc}"
      end

    end
  end
end