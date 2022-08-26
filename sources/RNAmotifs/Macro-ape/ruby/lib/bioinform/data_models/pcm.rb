require_relative '../support'
require_relative '../data_models'
require_relative '../conversion_algorithms/pcm2ppm_converter'
require_relative '../conversion_algorithms/pcm2pwm_converter'

module Bioinform
  class PCM < PM
    def count
      matrix.first.inject(&:+)
    end

    def to_pcm
      self
    end

    def to_pwm(pseudocount = Math.log(count))
      ConversionAlgorithms::PCM2PWMConverter.convert(self, pseudocount: pseudocount)
    end

    def to_ppm
      ConversionAlgorithms::PCM2PPMConverter.convert(self)
    end
  end
end