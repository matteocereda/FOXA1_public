class TransfacFormatter
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
      "ID #{name}\nBF StubSpeciesName\nP0\tA\tC\tG\tT\n"
    else
      raise 'Transfac should have the name field'
    end
  end
  
  def matrix_string
    motif.each_position.map.with_index{|pos,ind| 
      line_number = ind.to_s
      line_number = (line_number.size == 1) ? "0#{line_number}" : line_number
      line_number + ' ' + pos.join("\t")
    }.join("\n")
  end
  
  def footer
    #"XX\n//\n"
    "\nXX\n//"
  end
  
  def to_s
    header + matrix_string + footer
  end
end