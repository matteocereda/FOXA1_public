module Bioinform
  module ConversionAlgorithms
    module PCM2PPMConverter
    
      # parameters hash is ignored
      def self.convert(pcm, parameters = {})
        matrix = pcm.each_position.map do |pos|
          pos.map do |el|
            el.to_f / pcm.count
          end
        end
        PPM.new(pcm.get_parameters.merge(matrix: matrix))
      end
    end
  end
end


      