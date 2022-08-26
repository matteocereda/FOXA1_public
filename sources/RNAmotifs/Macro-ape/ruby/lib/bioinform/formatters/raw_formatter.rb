class RawFormatter
  attr_accessor :motif, :options
  
  def initialize(motif, options = {})
    @motif = motif
    
    default_options = {with_name: true, letters_as_rows: false}
    @options = default_options.merge(options)
  end
  
  def name
    motif.name
  end
  
  def header
    if options[:with_name] && name
      name + "\n"
    else
      ''
    end
  end
  
  def matrix_string
    if options[:letters_as_rows]
      hsh = motif.to_hash
      [:A,:C,:G,:T].collect{|letter| "#{letter}|" + hsh[letter].join("\t")}.join("\n")
    else
      motif.each_position.map{|pos| pos.join("\t")}.join("\n")
    end
  end
  
  def footer
    # "\n"
    ''
  end
  
  
  def to_s
    header + matrix_string + footer
  end
end