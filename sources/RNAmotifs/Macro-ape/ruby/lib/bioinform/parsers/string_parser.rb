require 'strscan'
require_relative '../support'
require_relative '../parsers/parser'

module Bioinform
  class StringParser < Parser
    include MultipleMotifsParser
    attr_reader :scanner, :row_acgt_markers

    def initialize(input)
      raise ArgumentError, 'StringParser should be initialized with a String'  unless input.is_a?(String)
      super
      @scanner = StringScanner.new(input.multiline_squish)
    end

    def number_pat
      /[+-]?\d+(\.\d+)?([eE][+-]?\d{1,3})?/
    end

    def header_pat
      />?\s*(?<name>\S+)\n/
    end

    def row_pat
      /([ACGT]\s*[:|]?\s*)?(?<row>(#{number_pat} )*#{number_pat})\n?/
    end

    def scan_row
      match = scanner.advanced_scan(row_pat)
      match && match[:row]
    end

    def split_row(row_string)
      row_string.split.map(&:to_f)
    end

    def scan_any_spaces
      scanner.scan(/\s+/)
    end

    def parse_name
      match = scanner.advanced_scan(header_pat)
      match && match[:name]
    end

    def parse_matrix
      matrix = []
      @row_acgt_markers = true  if scanner.check(/A.*\nC.*\nG.*\nT.*\n?/)
      while row_string = scan_row
        matrix << split_row(row_string)
      end
      matrix
    end

    def parse_acgt_header
      scanner.scan(/A\s*C\s*G\s*T\s*\n/i)
    end

    def parse!
      scan_any_spaces
      name = parse_name
      parse_acgt_header
      matrix = parse_matrix
      matrix = matrix.transpose if row_acgt_markers
      Parser.parse!(matrix).tap{|result| result.name = name}
    end

    def scanner_reset
      scanner.reset
    end
  end
end