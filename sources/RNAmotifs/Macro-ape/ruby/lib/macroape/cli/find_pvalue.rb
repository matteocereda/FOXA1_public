require_relative '../../macroape'

module Macroape
  module CLI
    module FindPValue

      def self.main(argv)
        doc = <<-EOS.strip_doc
          Command-line format:
          #{run_tool_cmd} <pat-file> <threshold list>... [options]

          Options:
            [-d <discretization level>]
            [--pcm] - treat the input file as Position Count Matrix. PCM-to-PWM transformation to be done internally.
            [-b <background probabilities] ACGT - 4 numbers, comma-delimited(spaces not allowed), sum should be equal to 1, like 0.25,0.24,0.26,0.25

          Examples:
            #{run_tool_cmd} motifs/KLF4_f2.pat 7.32
            #{run_tool_cmd} motifs/KLF4_f2.pat 7.32 4.31 5.42 -d 1000 -b 0.2,0.3,0.3,0.2
        EOS

        if argv.empty? || ['-h', '--h', '-help', '--help'].any?{|help_option| argv.include?(help_option)}
          STDERR.puts doc
          exit
        end

        discretization = 10000
        background = [1,1,1,1]
        thresholds = []
        max_hash_size = 10000000

        data_model = argv.delete('--pcm') ? Bioinform::PCM : Bioinform::PWM
        filename = argv.shift

        loop do
          begin
            Float(argv.first)
            thresholds << argv.shift.to_f
          rescue
            raise StopIteration
          end
        end

        raise 'No input. You should specify input file' unless filename
        raise 'You should specify at least one threshold' if thresholds.empty?

        until argv.empty?
          case argv.shift
            when '-b'
              background = argv.shift.split(',').map(&:to_f)
            when '-d'
              discretization = argv.shift.to_f
            when '--max-hash-size'
              max_hash_size = argv.shift.to_i
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

        counts = pwm.counts_by_thresholds(* thresholds.map{|count| count * discretization})
        infos = []
        thresholds.each do |threshold|
          count = counts[threshold * discretization]
          pvalue = count.to_f / pwm.vocabulary_volume
          infos << {threshold: threshold,
                    number_of_recognized_words: count,
                    pvalue: pvalue}
        end

        puts Helper.find_pvalue_info_string( infos,
                                            {discretization: discretization,
                                            background: background} )
      rescue => err
        STDERR.puts "\n#{err}\n#{err.backtrace.first(5).join("\n")}\n\nUse --help option for help\n\n#{doc}"
      end

    end
  end
end
