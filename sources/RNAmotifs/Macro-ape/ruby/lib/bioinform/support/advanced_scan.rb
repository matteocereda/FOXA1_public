require 'strscan'

class StringScanner
  def advanced_scan(pat)
    result = scan(pat)
    result && result.match(pat)
  end
end