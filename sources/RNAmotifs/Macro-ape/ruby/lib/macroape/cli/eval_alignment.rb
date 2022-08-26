require_relative '../../macroape'

module Macroape
  module CLI
    module EvalAlignment

      def self.main(argv)
        doc = <<-EOS.strip_doc
          Command-line format:
          #{run_tool_cmd} <1st matrix pat-file> <2nd matrix pat-file> <shift> <orientation(direct/revcomp)> [options]

          Options:
            [-p <P-value>]
            [-d <discretization level>]
            [--pcm] - treat the input file as Position Count Matrix. PCM-to-PWM transformation to be done internally.
            [--boundary lower|upper] Upper boundary (default) means that the obtained P-value is greater than or equal to the requested P-value
            [-b <background probabilities] ACGT - 4 numbers, comma-delimited(spaces not allowed), sum should be equal to 1, like 0.25,0.24,0.26,0.25
            [--first-threshold <threshold for the first matrix>]
            [--second-threshold <threshold for the second matrix>]

          Examples:
            #{run_tool_cmd} motifs/KLF4_f2.pat motifs/SP1_f1.pat -1 direct -p 0.0005 -d 100 -b 0.3,0.2,0.2,0.3
            #{run_tool_cmd} motifs/KLF4_f2.pat motifs/SP1_f1.pat 3 revcomp
        EOS

        if argv.empty? || ['-h', '--h', '-help', '--help'].any?{|help_option| argv.include?(help_option)}
          STDERR.puts doc
          exit
        end

        pvalue = 0.0005
        discretization = 10.0

        first_background = [1,1,1,1]
        second_background = [1,1,1,1]
        max_hash_size = 10000000
        max_pair_hash_size = 10000
        pvalue_boundary = :upper

        data_model = argv.delete('--pcm') ? Bioinform::PCM : Bioinform::PWM

        first_file = argv.shift
        second_file = argv.shift

        shift = argv.shift
        orientation = argv.shift

        raise 'You should specify two input sources (each is filename or .stdin)'  unless first_file and second_file
        raise 'You should specify shift' unless shift
        raise 'You should specify orientation' unless orientation

        shift = shift.to_i
        orientation = orientation.to_sym

        case orientation
          when :direct
            reverse = false
          when :revcomp
            reverse = true
          else
            raise 'Unknown orientation(direct/revcomp)'
        end


        until argv.empty?
          case argv.shift
            when '-p'
              pvalue = argv.shift.to_f
            when '-d'
              discretization = argv.shift.to_f
            when '--max-hash-size'
              max_hash_size = argv.shift.to_i
            when '--max-2d-hash-size'
              max_pair_hash_size = argv.shift.to_i
            when '-b'
              second_background = first_background = argv.shift.split(',').map(&:to_f)
            when '-b1'
              first_background = argv.shift.split(',').map(&:to_f)
            when '-b2'
              second_background = argv.shift.split(',').map(&:to_f)
            when '--boundary'
              pvalue_boundary = argv.shift.to_sym
              raise 'boundary should be either lower or upper'  unless  pvalue_boundary == :lower || pvalue_boundary == :upper
            when '--first-threshold'
              predefined_threshold_first = argv.shift.to_f
            when '--second-threshold'
              predefined_threshold_second = argv.shift.to_f
          end
        end
        raise 'background should be symmetric: p(A)=p(T) and p(G) = p(C)' unless first_background == first_background.reverse
        raise 'background should be symmetric: p(A)=p(T) and p(G) = p(C)' unless second_background == second_background.reverse

        if first_file == '.stdin' || second_file == '.stdin'
          input = $stdin.read
          parser = data_model.choose_parser(input).new(input)
        end

        if first_file == '.stdin'
          input_first = parser.parse
        else
          raise "Error! File #{first_file} don't exist" unless File.exist?(first_file)
          input_first = File.read(first_file)
        end
        pwm_first = data_model.new(input_first).to_pwm

        if second_file == '.stdin'
          input_second = parser.parse
        else
          raise "Error! File #{second_file} don't exist" unless File.exist?(second_file)
          input_second = File.read(second_file)
        end
        pwm_second = data_model.new(input_second).to_pwm

        pwm_first.set_parameters(background: first_background, max_hash_size: max_hash_size).discrete!(discretization)
        pwm_second.set_parameters(background: second_background, max_hash_size: max_hash_size).discrete!(discretization)

        cmp = Macroape::PWMCompareAligned.new(pwm_first, pwm_second, shift, orientation).set_parameters(max_pair_hash_size: max_pair_hash_size)

        if predefined_threshold_first
          threshold_first = predefined_threshold_first * discretization
        else
          if pvalue_boundary == :lower
            threshold_first = pwm_first.threshold(pvalue)
          else
            threshold_first = pwm_first.weak_threshold(pvalue)
          end
        end

        if predefined_threshold_second
          threshold_second = predefined_threshold_second * discretization
        else
          if pvalue_boundary == :lower
            threshold_second = pwm_second.threshold(pvalue)
          else
            threshold_second = pwm_second.weak_threshold(pvalue)
          end
        end
        info = cmp.alignment_infos.merge( cmp.jaccard(threshold_first, threshold_second) )
        info.merge!(predefined_threshold_first: predefined_threshold_first,
                    predefined_threshold_second: predefined_threshold_second,
                    threshold_first: threshold_first / discretization,
                    threshold_second: threshold_second / discretization,
                    discretization: discretization,
                    first_background: first_background,
                    second_background: second_background,
                    requested_pvalue: pvalue,
                    pvalue_boundary: pvalue_boundary)
        puts Helper.similarity_info_string(info)

      rescue => err
        STDERR.puts "\n#{err}\n#{err.backtrace.first(5).join("\n")}\n\nUse --help option for help\n\n#{doc}"
      end

    end
  end
end