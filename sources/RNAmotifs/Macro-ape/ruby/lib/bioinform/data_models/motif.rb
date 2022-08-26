require 'ostruct'
require 'active_support/core_ext/object/try'
require_relative '../support/parameters'
module Bioinform
  class Motif
    include Parameters
    make_parameters :pcm, :pwm, :ppm, :name, :original_data_model

    # 0)Motif.new()
    # 1)Motif.new(pcm: ..., pwm: ..., name: ...,threshold: ...)
    # 2)Motif.new(my_pcm)
    # 3)Motif.new(pm: my_pcm, threshold: ...)
    # 2) and 3) cases will automatically choose data model
    #### What if pm already is a Motif
    def initialize(parameters = {})
      case parameters
      when PM
        pm = parameters
        motif_type = pm.class.name.downcase.sub(/^.+::/,'').to_sym
        self.original_data_model = motif_type
        set_parameters(motif_type => pm)
      when Hash
        if parameters.has_key?(:pm) && parameters[:pm].is_a?(PM)
          pm = parameters.delete(:pm)
          motif_type = pm.class.name.downcase.sub(/^.+::/,'').to_sym
          self.original_data_model = motif_type
          set_parameters(motif_type => pm)
        end
        set_parameters(parameters)
      else
        raise ArgumentError, "Motif::new doesn't accept argument #{parameters} of class #{parameters.class}"
      end
    end

    def pm; ((original_data_model || :pm) == :pm) ? parameters.pm : send(original_data_model); end
    #def pcm; parameters.pcm; end
    def pwm; parameters.pwm || pcm.try(:to_pwm); end
    def ppm; parameters.ppm || pcm.try(:to_ppm); end
    #def pcm=(pcm); parameters.pcm = pcm; end
    #def pwm=(pwm); parameters.pwm = pwm; end
    #def ppm=(ppm); parameters.ppm = ppm; end
    def name; parameters.name || pm.name; end

    def method_missing(meth, *args)
      parameters.__send__(meth, *args)
    end

    def ==(other)
      parameters == other.parameters
    end
    
    def to_s
      parameters.to_s
    end
  end
end