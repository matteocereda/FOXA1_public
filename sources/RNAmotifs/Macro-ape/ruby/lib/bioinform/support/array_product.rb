class Array
  def self.product(*arrays)
    return []  if arrays.empty?
    arrays.first.product(*arrays[1..-1])
  end
end