require 'ostruct'
require_relative 'motif'

module Bioinform
  class Collection
    attr_accessor :container

    include Parameters
    make_parameters :name

    # collection name is a tag name for each motif in a collection. But motif can be included in several collections so have several tags
    def initialize(parameters = {})
      @container = []
      @parameters = OpenStruct.new(parameters)
      yield @parameters  if block_given?
    end

    def size
      container.size
    end

    def to_s(with_name = true)
      result = (with_name) ? "Collection: #{name.to_s}\n" : ''
      each do |pm, infos|
        result << pm.to_s << "\n\n"
      end
      result
    end

    def +(other)
      result = self.class.new
      container.each do |motif|
        result.container << motif
      end
      other.container.each do |motif|
        result.container << motif
      end
      result
    end

    def add_pm(pm, info)
#      pm.mark(self)
      container << Motif.new(info.marshal_dump.merge(pm: pm))
      #### What if pm already is a Motif
      self
    end

    def <<(pm)
      add_pm(pm, OpenStruct.new)
    end

    # collection.each{|motif| ... }
    # collection.each(:pwm, :threshold){|pwm,threshold| }
    def each(*args)
      if block_given?
        if args.empty?
          container.each{|motif| yield motif}
        else
          container.each{|motif| yield( *args.map{|arg| motif.parameters.send(arg)} ) }
        end
      else
        Enumerator.new(self, :each, *args)
      end
    end

    include Enumerable

    def ==(other)
      (parameters == other.parameters) && (container == other.container)
    rescue
      false
    end

  end
end