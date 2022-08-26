require 'ostruct'
require_relative '../support'
require_relative '../data_models/pm'
require_relative 'splittable_parser'

module Bioinform
  class Error < StandardError; end
  class ParsingError < Error; end
  class InvalidMatrix < Error; end
  
  class Parser
    attr_reader :input

    def initialize(*input)
      if input.size == 1  # [ [1,2,3,4] ],  [  [[1,2,3,4],[5,6,7,8]] ]
        if input.first.is_a?(Array) && input.first.all?{|el| el.is_a? Numeric}  # [ [1,2,3,4] ]
          @input = input
        else  # [  [[1,2,3,4],[5,6,7,8]] ]
          @input = input.first
        end
      else #[ [1,2,3,4], [5,6,7,8] ], [   ]
        @input = input
      end
    end

    def parse!
      matrix = self.class.transform_input(input)
      raise InvalidMatrix unless self.class.valid_matrix?(matrix)
      OpenStruct.new(matrix: matrix)
    end

    def parse
      parse! rescue nil
    end

    def self.choose(input, data_model = PM)
      data_model.choose_parser(input).new(input)
    end

    def self.parse!(*input)
      self.new(*input).parse!
    end
    def self.parse(*input)
      self.new(*input).parse
    end

    def self.valid_matrix?(matrix)
      PM.valid_matrix?(matrix)
    end

    # {A: 1, C: 2, G: 3, T: 4}  -->  [1,2,3,4]
    # {A: [1,2], C: [3,4], G: [5,6], T: [7,8]}  --> [[1,3,5,7],[2,4,6,8]] ( == [[1,2], [3,4], [5,6], [7,8]].transpose)
    def self.array_from_acgt_hash(hsh)
      hsh = normalize_hash_keys(hsh)
      raise 'some of hash keys A,C,G,T are missing or hash has excess keys' unless hsh.keys.sort == [:A,:C,:G,:T]
      result = [:A,:C,:G,:T].collect{|letter| hsh[letter] }
      result.all?{|el| el.is_a?(Array)} ? result.transpose : result
    end

    # {a: 1, C: 2, 'g' => 3, 'T' => 4} --> {A: 1, C: 2, G: 3, T: 4}
    def self.normalize_hash_keys(hsh)
      hsh.collect_hash{|key,value| [key.to_s.upcase.to_sym, value] }
    end

    # [[1,2,3,4], [2,3,4,5]] --> [[1,2,3,4], [2,3,4,5]]
    # [{A:1, C:2, G:3, T:4}, {A:2, C:3, G:4, T:5}] --> [{A:1, C:2, G:3, T:4}, {A:2, C:3, G:4, T:5}]
    # {:A => [1,2,3], :c => [2,3,4], 'g' => [3,4,5], 'T' => [4,5,6]} --> [[1,2,3],[2,3,4],[3,4,5],[4,5,6]].transpose
    def self.try_convert_to_array(input)
      case input
      when Array then input
      when Hash then array_from_acgt_hash(input)
      else raise TypeError, 'input of Bioinform::Parser::array_from_acgt_hash should be Array or Hash'
      end
    end

    def self.transform_input(input)
      result = try_convert_to_array(input).map{|el| try_convert_to_array(el)}
      need_tranpose?(result) ? result.transpose : result
    end

    # point whether matrix input positions(need not be transposed -- false) or letters(need -- true) as first index
    # [[1,3,5,7], [2,4,6,8]] --> false
    # [[1,2],[3,4],[5,6],[7,8]] --> true
    def self.need_tranpose?(input)
      (input.size == 4) && input.any?{|x| x.size != 4}
    end
  end
end