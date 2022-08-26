module Enumerable
  # %w{A C G T}.collect_hash{|k| [k*2, k*3] }
  # # ==> {"AA" => "AAA", "CC" => "CCC", "GG" => "GGG", "TT" => "TTT"}
  def collect_hash(&block)
    block_given?  ?  Hash[ collect(&block) ]  :  Hash[ collect{|k,v| [k,v]} ]
  end
end