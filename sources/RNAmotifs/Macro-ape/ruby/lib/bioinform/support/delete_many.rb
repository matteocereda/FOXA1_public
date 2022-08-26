class Array
  def delete_at_many(*indices)
    indices.uniq.sort.reverse.each{|ind| delete_at ind}
  end
  def delete_many(*elements)
    elements.each{|el| delete el}
  end
end

class Hash
  def delete_many(*keys)
    keys.each{|el| delete el}
  end
end