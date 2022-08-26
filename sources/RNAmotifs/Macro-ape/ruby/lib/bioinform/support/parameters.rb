require 'ostruct'
module Bioinform
  module Parameters
    def self.included(base)
      base.extend(ClassMethods)
    end
    module ClassMethods
      def make_parameters(*params)
        params.each do |param|
          define_method(param){ parameters.send(param) }
          define_method("#{param}="){|new_value| parameters.send("#{param}=", new_value) }
        end
      end
    end
    def parameters; @parameters ||= OpenStruct.new; end
    def set_parameters(hsh)
      hsh.each{|k,v| send("#{k}=", v) }
      self
    end
    # return hash of parameters
    def get_parameters
      @parameters.marshal_dump
    end
    def parameter_defined?(param_name)
      get_parameters.has_key?(param_name)
    end
  end
end