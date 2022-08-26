require_relative '../support'
require_relative '../data_models'
module Bioinform
  class PWM < PM
    def score_mean
      each_position.inject(0){ |mean, position| mean + position.each_index.inject(0){|sum, letter| sum + position[letter] * probability[letter]} }
    end
    def score_variance
      each_position.inject(0) do |variance, position|
        variance  + position.each_index.inject(0) { |sum,letter| sum + position[letter]**2 * probability[letter] } -
                    position.each_index.inject(0) { |sum,letter| sum + position[letter]    * probability[letter] }**2
      end
    end

    def threshold_gauss_estimation(pvalue)
      sigma = Math.sqrt(score_variance)
      n_ = Math.inverf(1 - 2 * pvalue) * Math.sqrt(2)
      score_mean + n_ * sigma
    end

    def score(word)
      word = word.upcase
      raise ArgumentError, 'word in PWM#score(word) should have the same length as matrix'  unless word.length == length
      #raise ArgumentError, 'word in PWM#score(word) should have only ACGT-letters'  unless word.each_char.all?{|letter| %w{A C G T}.include? letter}
      (0...length).map do |pos|
        begin
        # Need support of N-letters and other IUPAC
          letter = word[pos]
          matrix[pos][IndexByLetter[letter]]
        rescue
          raise ArgumentError, 'word in PWM#score(word) should have only ACGT-letters'
        end
      end.inject(&:+)
    end

    def to_pwm
      self
    end

    def best_score
      @matrix.inject(0.0){|sum, col| sum + col.max}
    end
    def worst_score
      @matrix.inject(0.0){|sum, col| sum + col.min}
    end

    # best score of suffix s[i..l]
    def best_suffix(i)
      @matrix[i...length].map(&:max).inject(0.0, &:+)
    end

    def worst_suffix(i)
      @matrix[i...length].map(&:min).inject(0.0, &:+)
    end
  end
end