require 'ostruct'
require_relative '../support'
require_relative '../parsers'
require_relative '../formatters'

module Bioinform
  IndexByLetter = {'A' => 0, 'C' => 1, 'G' => 2, 'T' => 3, A: 0, C: 1, G: 2, T: 3}
  LetterByIndex = {0 => :A, 1 => :C, 2 => :G, 3 => :T}

  class PM
    attr_accessor :matrix, :parameters

    include Parameters
    make_parameters  :name, :background # , :tags

#    def mark(tag)
#      tags << tag
#    end

#    def tagged?(tag)
#      tags.any?{|t| (t.eql? tag) || (t.respond_to?(:name) && t.name && (t.name == tag)) }
#    end

    def self.choose_parser(input)
      [TrivialParser, YAMLParser, Parser, StringParser, StringFantomParser, JasparParser, TrivialCollectionParser, YAMLCollectionParser].find do |parser|
        self.new(input, parser) rescue nil
      end
    end

    def self.split_on_motifs(input)
      parser = choose_parser(input)
      raise ParsingError, "No parser can parse given input"  unless parser
      parser.split_on_motifs(input, self)
    end

    def initialize(input, parser = nil)
      @parameters = OpenStruct.new
      parser ||= self.class.choose_parser(input)
      raise 'No one parser can process input'  unless parser
      result = parser.new(input).parse
      @matrix = result.matrix
      self.name = result.name
#      self.tags = result.tags || []
      self.background = result.background || [1, 1, 1, 1]
      raise 'matrix not valid'  unless valid?
    end

    def ==(other)
      @matrix == other.matrix && background == other.background && name == other.name
    rescue
      false
    end

    def self.valid_matrix?(matrix)
      matrix.is_a?(Array) &&
      ! matrix.empty? &&
      matrix.all?{|pos| pos.is_a?(Array)} &&
      matrix.all?{|pos| pos.size == 4} &&
      matrix.all?{|pos| pos.all?{|el| el.is_a?(Numeric)}}
    rescue
      false
    end

    def valid?
      self.class.valid_matrix?(@matrix)
    end

    def each_position
      if block_given?
        matrix.each{|pos| yield pos}
      else
        Enumerator.new(self, :each_position)
      end
    end

    def length
      @matrix.length
    end
    alias_method :size, :length

    def to_s(options = {}, formatter = RawFormatter)
      formatter.new(self, options).to_s
    end

    def pretty_string(options = {})
      default_options = {with_name: true, letters_as_rows: false}

      return to_s(options)  if options[:letters_as_rows]

      options = default_options.merge(options)
      header = %w{A C G T}.map{|el| el.rjust(4).ljust(7)}.join + "\n"
      matrix_rows = each_position.map do |position|
        position.map{|el| el.round(3).to_s.rjust(6)}.join(' ')
      end

      matrix_str = matrix_rows.join("\n")

      if options[:with_name] && name
        name + "\n" + header + matrix_str
      else
        header + matrix_str
      end
    end

    def to_hash
      hsh = %w{A C G T}.each_with_index.collect_hash do |letter, letter_index|
        [ letter, @matrix.map{|pos| pos[letter_index]} ]
      end
      hsh.with_indifferent_access
    end

    def self.zero_column
      [0, 0, 0, 0]
    end

    def reverse_complement!
      @matrix.reverse!.map!(&:reverse!)
      self
    end
    def left_augment!(n)
      n.times{ @matrix.unshift(self.class.zero_column) }
      self
    end
    def right_augment!(n)
      n.times{ @matrix.push(self.class.zero_column) }
      self
    end

    def discrete!(rate)
      @matrix.map!{|position| position.map{|element| (element * rate).ceil}}
      self
    end

    def vocabulary_volume
      background.inject(&:+) ** length
    end

    def probability
      sum = background.inject(0.0, &:+)
      background.map{|element| element.to_f / sum}
    end

    def reverse_complement
      dup.reverse_complement!
    end
    def left_augment(n)
      dup.left_augment!(n)
    end
    def right_augment(n)
      dup.right_augment!(n)
    end
    def discrete(rate)
      dup.discrete!(rate)
    end
    def dup
      deep_dup
    end

    def as_pcm
      PCM.new(get_parameters.merge(matrix: matrix))
    end
    def as_ppm
      PPM.new(get_parameters.merge(matrix: matrix))
    end
    def as_pwm
      PWM.new(get_parameters.merge(matrix: matrix))
    end
  end
end