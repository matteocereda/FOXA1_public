require_relative '../../macroape'

module Macroape
  module CLI
    module FindThreshold

      def self.main(argv)
        doc = <<-EOS.strip_doc
          Command-line format:
          #{run_tool_cmd} <pat-file> [<list of P-values>...] [options]

          Options:
            [-d <discretization level>]
            [--pcm] - treat the input file as Position Count Matrix. PCM-to-PWM transformation to be done internally.
            [--boundary lower|upper] Lower boundary (default) means that the obtained P-value is less than or equal to the requested P-value
            [-b <background probabilities] ACGT - 4 numbers, comma-delimited(spaces not allowed), sum should be equal to 1, like 0.25,0.24,0.26,0.25

          Example:
            #{run_tool_cmd} motifs/KLF4_f2.pat
            #{run_tool_cmd} motifs/KLF4_f2.pat 0.001 0.0001 0.0005 -d 1000 -b 0.4,0.3,0.2,0.1
        EOS

        if argv.empty? || ['-h', '--h', '-help', '--help'].any?{|help_option| argv.include?(help_option)}
          STDERR.puts doc
          exit
        end

        background = [1,1,1,1]
        default_pvalues = [0.0005]
        discretization = 10000
        max_hash_size = 10000000
        data_model = argv.delete('--pcm') ? Bioinform::PCM : Bioinform::PWM

        pvalue_boundary = :lower


        filename = argv.shift
        raise 'No input. You should specify input file' unless filename

        pvalues = []
        loop do
          begin
            Float(argv.first)
            pvalues << argv.shift.to_f
          rescue
            raise StopIteration
          end
        end
        pvalues = default_pvalues if pvalues.empty?

        until argv.empty?
          case argv.shift
            when '-b'
              background = argv.shift.split(',').map(&:to_f)
            when '--max-hash-size'
              max_hash_size = argv.shift.to_i
            when '-d'
              discretization = argv.shift.to_f
            when '--boundary'
              pvalue_boundary = argv.shift.to_sym
              raise 'boundary should be either lower or upper'  unless  pvalue_boundary == :lower || pvalue_boundary == :upper
            end
        end

        if filename == '.stdin'
          input = $stdin.read
        else
          raise "Error! File #{filename} doesn't exist" unless File.exist?(filename)
          input = File.read(filename)
        end
        pwm = data_model.new(input).to_pwm
        pwm.set_parameters(background: background, max_hash_size: max_hash_size).discrete!(discretization)

        infos = []
        collect_infos_proc = ->(pvalue, threshold, real_pvalue) do
          infos << {expected_pvalue: pvalue,
                    threshold: threshold / discretization,
                    real_pvalue: real_pvalue,
                    recognized_words: pwm.vocabulary_volume * real_pvalue }
        end
        if pvalue_boundary == :lower
          pwm.thresholds(*pvalues, &collect_infos_proc)
        else
          pwm.weak_thresholds(*pvalues, &collect_infos_proc)
        end
        puts Helper.threshold_infos_string(infos,
                                          {discretization: discretization,
                                          background: background,
                                          pvalue_boundary: pvalue_boundary} )
      rescue => err
        STDERR.puts "\n#{err}\n#{err.backtrace.first(5).join("\n")}\n\nUse --help option for help\n\n#{doc}"
      end

    end
  end
end