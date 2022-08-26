require_relative '../support'
require_relative 'parser'
require_relative '../data_models/collection'
require 'yaml'

module Bioinform
  # YAMLParser can be used to parse hashes returned by #parse method of other parsers:
  # yaml_dump_of_pm = PM.new(...).to_yaml
  # PM.new(yaml_dump_of_pm, YAMLParser)
  class YAMLParser < Parser
    def initialize(input)
      @input = input
    end
    def parse!
      YAML.load(input)
    rescue Psych::SyntaxError
      raise 'parsing error'
    end
  end

  class YAMLCollectionParser < Parser
    include MultipleMotifsParser
    def initialize(input)
      @input = input
    end
    def collection
      @collection ||= YAML.load(input)
    end
    def parse!
      collection.container.shift.pm
    rescue Psych::SyntaxError
      raise 'parsing error'
    end
  end
end