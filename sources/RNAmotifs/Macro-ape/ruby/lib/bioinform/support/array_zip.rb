class Array
  def self.zip(*arrays)
    return []  if arrays.empty?
    arrays.first.zip(*arrays[1..-1])
  end
end