module Macroape
  class PWMCompareAligned
    # unoptimized version of this and related methods
    def counts_for_two_matrices(threshold_first, threshold_second)
      # just not to call method each time
      first_background = first.background
      second_background = second.background
      unless first_background == second_background
        first_result = get_counts(threshold_first, threshold_second) {|score,letter| first_background[letter] * score }
        second_result = get_counts(threshold_first, threshold_second) {|score,letter| second_background[letter] * score }
        return [first_result, second_result]
      end
      if first.background == [1,1,1,1]
        result = get_counts(threshold_first, threshold_second) {|score,letter| score}
        [result, result]
      else
        result = get_counts(threshold_first, threshold_second) {|score,letter| first_background[letter] * score }
        [result, result]
      end
    end


    # block has form: {|score,letter| contribution to count by `letter` with `score` }
    def get_counts(threshold_first, threshold_second, &count_contribution_block)
      # scores_on_first_pwm, scores_on_second_pwm --> count
      scores = { 0 => {0 => 1} }
      length.times do |column|
        new_scores = recalc_score_hash(scores,
                          first.matrix[column], second.matrix[column],
                          threshold_first - first.best_suffix(column + 1),
                          threshold_second - second.best_suffix(column + 1), &count_contribution_block)
        scores.replace(new_scores)
        if max_pair_hash_size && scores.inject(0){|sum,hsh|sum + hsh.size} > max_pair_hash_size
          raise 'Hash overflow in Macroape::AlignedPairIntersection#counts_for_two_matrices_with_different_probabilities'
        end
      end
      scores.inject(0.0){|sum,(score_first, hsh)| sum + hsh.inject(0.0){|sum,(score_second, count)| sum + count }}
    end

    # wouldn't work without count_contribution_block
    def recalc_score_hash(scores, first_column, second_column, least_sufficient_first, least_sufficient_second)
      new_scores = Hash.new{|h,k| h[k] = Hash.new(0)}
      scores.each do |score_first, second_scores|
        second_scores.each do |score_second, count|

          4.times do |letter|
            new_score_first = score_first + first_column[letter]
            if new_score_first >= least_sufficient_first
              new_score_second = score_second + second_column[letter]
              if new_score_second >= least_sufficient_second
                new_scores[new_score_first][new_score_second] += yield(count, letter)
              end
            end
          end

        end
      end
      new_scores
    end

  end
end