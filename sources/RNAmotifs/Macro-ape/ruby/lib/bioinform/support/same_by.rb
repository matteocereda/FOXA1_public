module Enumerable
  def same_by?(&block)
    return true if empty?
    if block_given?
      first_result = yield(first)
      all?{|el| first_result == yield(el)}
    else
      first_result = first
      all?{|el| first_result == el}
    end
  end
end
