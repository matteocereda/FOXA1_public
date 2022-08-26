module Bioinform
  class Parser
    module SingleMotifParser
      def self.included(base)
        base.class_eval { extend ClassMethods }
        include Enumerable
        alias_method :split, :to_a
      end
      module ClassMethods
        def split_on_motifs(input, pm_klass = PM)
          [ input.is_a?(pm_klass) ? self : pm_klass.new(input, self) ]
        end
      end
      def each
        if block_given?
          yield self
        else
          Enumerator.new(self, :each)
        end
      end
    end
    include SingleMotifParser

    module MultipleMotifsParser
      def self.included(base)
        base.class_eval { extend ClassMethods }
        include Enumerable
        alias_method :split, :to_a
      end
      module ClassMethods
        def split_on_motifs(input, pm_klass = PM)
          split(input).map{|el| el.is_a?(pm_klass) ? el : pm_klass.new(el)}
        end
        def split(input)
          self.new(input).split
        end
        private :split
      end

      def scanner_reset
      end

      def each
        if block_given?
          scanner_reset
          while result = parse
            yield result
          end
        else
          Enumerator.new(self, :each)
        end
      end

      private :scanner_reset
    end
  end
end