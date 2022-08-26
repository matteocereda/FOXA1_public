def strip_doc(doc)
  doc.strip_doc
end

class String
  def strip_doc
    gsub(/^#{self[/\A +/]}/,'')
  end
end