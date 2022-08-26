require_relative 'collect_hash'

class Array
  def partial_sums(initial = 0.0)
    sums = initial
    map{|el| sums += el}
  end
end

class Hash
# {1 => 5, 4 => 3, 3 => 2}.partial_sums == {1=>5, 3=>7, 4=>10}
  def partial_sums(initial = 0.0)
    sums = initial
    sort.collect_hash{|k,v| [k, sums += v]}
  end
end