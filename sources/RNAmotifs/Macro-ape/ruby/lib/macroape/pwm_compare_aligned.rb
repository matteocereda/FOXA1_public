require 'bioinform/support/parameters'
require_relative 'aligned_pair_intersection'

module Macroape
  class PWMCompareAligned
    include Bioinform::Parameters
    # sets or gets limit of summary size of calculation hash. It's a defence against overuse CPU resources by non-appropriate data
    make_parameters :max_pair_hash_size

    attr_reader :first, :second, :length, :shift, :orientation, :first_length, :second_length, :parameters

    def initialize(first_unaligned, second_unaligned, shift, orientation)
      @parameters = OpenStruct.new
      @shift, @orientation = shift, orientation

      @first_length, @second_length = first_unaligned.length, second_unaligned.length
      @length = self.class.calculate_alignment_length(@first_length, @second_length, @shift)

      first, second = first_unaligned, second_unaligned

      if shift > 0
        second = second.left_augment(shift)
      else
        first = first.left_augment(-shift)
      end

      @first = first.right_augment(@length - first.length)
      @second = second.right_augment(@length - second.length)
    end

    def direct?
      orientation == :direct
    end

    def overlap
      length.times.count{|pos| first_overlaps?(pos) && second_overlaps?(pos) }
    end

    def first_pwm_alignment
      length.times.map do |pos|
        if first_overlaps?(pos)
          '>'
        else
          '.'
        end
      end.join
    end

    def second_pwm_alignment
      length.times.map do |pos|
        if second_overlaps?(pos)
          direct? ? '>' : '<'
        else
          '.'
        end
      end.join
    end

    def alignment_infos
      {shift: shift,
      orientation: orientation,
      text: "#{first_pwm_alignment}\n#{second_pwm_alignment}",
      overlap: overlap,
      alignment_length: length}
    end

    # whether first matrix overlap specified position of alignment
    def first_overlaps?(pos)
      return false unless pos >= 0 && pos < length
      if shift > 0
        pos < first_length
      else
        pos >= -shift && pos < -shift + first_length
      end
    end

    def second_overlaps?(pos)
      return false unless pos >= 0 && pos < length
      if shift > 0
        pos >= shift && pos < shift + second_length
      else
        pos < second_length
      end
    end

    def jaccard(first_threshold, second_threshold)
      f = first.count_by_threshold(first_threshold)
      s = second.count_by_threshold(second_threshold)
      if f == 0 || s == 0
        return {similarity: -1, tanimoto: -1, recognized_by_both: 0,
              recognized_by_first: f,
              recognized_by_second: s,
            }
      end

      intersect = counts_for_two_matrices(first_threshold, second_threshold)
      intersect = Math.sqrt(intersect[0] * intersect[1])
      union = f + s - intersect
      similarity = intersect.to_f / union
      { similarity: similarity,  tanimoto: 1.0 - similarity,  recognized_by_both: intersect,
        recognized_by_first: f,  recognized_by_second: s,
        real_pvalue_first: f / first.vocabulary_volume, real_pvalue_second: s / second.vocabulary_volume }
    end

    def jaccard_by_pvalue(pvalue)
      threshold_first = first.threshold(pvalue)
      threshold_second = second.threshold(pvalue)
      jaccard(threshold_first, threshold_second)
    end

    def jaccard_by_weak_pvalue(pvalue)
      threshold_first = first.weak_threshold(pvalue)
      threshold_second = second.weak_threshold(pvalue)
      jaccard(threshold_first, threshold_second)
    end

    def self.calculate_alignment_length(first_len, second_len, shift)
      if shift > 0
        [first_len, second_len + shift].max
      else
        [first_len - shift, second_len].max
      end
    end
  end

end