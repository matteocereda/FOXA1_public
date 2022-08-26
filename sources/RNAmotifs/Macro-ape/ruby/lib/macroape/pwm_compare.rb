require 'bioinform/support/parameters'

module Macroape
  class PWMCompare
    include Bioinform::Parameters
    # sets or gets limit of summary size of calculation hash. It's a defence against overuse CPU resources by non-appropriate data
    make_parameters :max_pair_hash_size

    attr_reader :first, :second, :parameters
    def initialize(first, second)
      @parameters = OpenStruct.new
      @first = first
      @second = second
    end

    def jaccard(threshold_first, threshold_second)
      self.map_each_alignment do |alignment|
        alignment.alignment_infos.merge( alignment.jaccard(threshold_first, threshold_second) )
      end.max_by {|alignment_infos| alignment_infos[:similarity] }
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

    def each_alignment
      (-second.length..first.length).to_a.product([:direct]) do |shift, orientation|
        yield PWMCompareAligned.new(first, second, shift, orientation).set_parameters(max_pair_hash_size: max_pair_hash_size)
      end
    end

    include Enumerable
    alias_method :each, :each_alignment
    alias_method :map_each_alignment, :map
  end
end