require 'bioinform/support/strip_doc'

class String
  def snake_case
    gsub(/[A-Z]+/){|big| "_#{big.downcase}" }.sub(/^_/,'')
  end
end

class Module
  def run_tool_cmd
    if Macroape::STANDALONE
      "ruby #{tool_name}.rb"
    else
      tool_name
    end
  end
  def tool_name
    self.name.split('::').last.snake_case
  end
end

module Macroape
  module CLI
    class OutputInformation
      def initialize(data = nil)
        @table_parameter_descriptions = []

        @parameter_descriptions = []
        @parameter_value_infos = []

        @resulting_value_descriptions = []
        @resulting_value_infos = []

        @table_headers = []
        @table_rows = []
        @table_rows_callbacks = []
        @data = data
        yield self  if block_given?
      end

      def parameters_info
        [*@parameter_descriptions, *@parameter_value_infos]
      end
      def resulting_values_info
        [*@resulting_value_descriptions, *@resulting_value_infos]
      end
      def result
        [parameters_info, resulting_values_info, resulting_table].reject(&:empty?).map{|b|b.join("\n")}.join("\n#\n")
        #[*parameters_info, '#', *resulting_values_info, '#', *resulting_table].join("\n")
      end

      def add_parameter(param_name, description, value, &block)
        @parameter_descriptions << parameter_description_string(param_name, description)
        @parameter_value_infos << "# #{param_name} = #{value}"
      end

      def add_resulting_value(param_name, description, value, &block)
        @resulting_value_descriptions << parameter_description_string(param_name, description)
        @resulting_value_infos << "#{param_name}\t#{value}"
      end

      def add_table_parameter(param_name, description, key_in_hash, &block)
        @table_parameter_descriptions << parameter_description_string(param_name, description)
        add_table_parameter_without_description(param_name, key_in_hash, &block)
      end

      def add_table_parameter_without_description(param_name, key_in_hash, &block)
        @table_headers << param_name
        @table_rows << key_in_hash
        @table_rows_callbacks << block
      end

      def parameter_description_string(param_name, description)
        "# #{param_name}: #{description}"
      end

      def table_content
        @data.map{|info|
          @table_rows.zip(@table_rows_callbacks).map{|row,callback| callback ? callback.call(info[row]) : info[row] }.join("\t")
        }
      end

      def header_content
        '# ' + @table_headers.join("\t")
      end

      def resulting_table
        @data ? [*@table_parameter_descriptions, header_content, *table_content] : []
      end

      # printed only if it is not wordwise [1,1,1,1]
      def background_parameter(param_name, description, value, &block)
        add_parameter(param_name, description, value.join(','), &block)  unless value == [1,1,1,1]
      end
    end

    module Helper

      def self.similarity_info_string(info)
        OutputInformation.new { |infos|
          infos.add_parameter('V', 'discretization', info[:discretization] )
          infos.add_parameter('P', 'requested P-value', info[:requested_pvalue])  unless info[:predefined_threshold_first] && info[:predefined_threshold_second]
          infos.add_parameter('T1', 'threshold for the 1st matrix', info[:predefined_threshold_first] )  if info[:predefined_threshold_first]
          infos.add_parameter('T2', 'threshold for the 2nd matrix', info[:predefined_threshold_second] )  if info[:predefined_threshold_second]
          infos.add_parameter('PB', 'P-value boundary', info[:pvalue_boundary])
          if info[:first_background] == info[:second_background]
            infos.background_parameter('B', 'background', info[:first_background])
          else
            infos.background_parameter('B1', 'background for the 1st model', info[:first_background])
            infos.background_parameter('B2', 'background for the 2nd model', info[:second_background])
          end

          infos.add_resulting_value('S', 'similarity', info[:similarity])
          infos.add_resulting_value('D', 'distance (1-similarity)', info[:tanimoto])
          infos.add_resulting_value('L', 'length of the alignment', info[:alignment_length])
          infos.add_resulting_value('SH', 'shift of the 2nd PWM relative to the 1st', info[:shift])
          infos.add_resulting_value('OR', 'orientation of the 2nd PWM relative to the 1st', info[:orientation])
          infos.add_resulting_value('A1', 'aligned 1st matrix', info[:text].lines.to_a.first.strip )
          infos.add_resulting_value('A2', 'aligned 2nd matrix', info[:text].lines.to_a.last.strip )
          infos.add_resulting_value('W', 'number of words recognized by both models (model = PWM + threshold)', info[:recognized_by_both] )
          infos.add_resulting_value('W1', 'number of words and recognized by the first model', info[:recognized_by_first] )
          infos.add_resulting_value('P1', 'P-value for the 1st matrix', info[:real_pvalue_first] )
          infos.add_resulting_value('T1', 'threshold for the 1st matrix', info[:threshold_first] )  unless info[:predefined_threshold_first]
          infos.add_resulting_value('W2', 'number of words recognized by the 2nd model', info[:recognized_by_second] )
          infos.add_resulting_value('P2', 'P-value for the 2nd matrix', info[:real_pvalue_second] )
          infos.add_resulting_value('T2', 'threshold for the 2nd matrix', info[:threshold_second] )  unless info[:predefined_threshold_second]
        }.result
      end

############################################

      def self.threshold_infos_string(data, parameters)
        OutputInformation.new(data) { |infos|
          infos.add_parameter('V', 'discretization value', parameters[:discretization])
          infos.add_parameter('PB', 'P-value boundary', parameters[:pvalue_boundary])
          infos.background_parameter('B', 'background', parameters[:background])

          infos.add_table_parameter('P', 'requested P-value', :expected_pvalue)
          infos.add_table_parameter('AP', 'actual P-value', :real_pvalue)
          infos.add_table_parameter('W', 'number of recognized words', :recognized_words)  if parameters[:background] == [1, 1, 1, 1]
          infos.add_table_parameter('T', 'threshold', :threshold)
        }.result
      end

############################################

      def self.scan_collection_infos_string(data, parameters)
        OutputInformation.new(data) { |infos|
          infos.add_parameter('MS', 'minimal similarity to output', parameters[:cutoff])
          infos.add_parameter('P', 'P-value', parameters[:pvalue])
          infos.add_parameter('PB', 'P-value boundary', parameters[:pvalue_boundary])
          if parameters[:precision_mode] == :precise
            infos.add_parameter('VR', 'discretization value, rough', parameters[:rough_discretization])
            infos.add_parameter('VP', 'discretization value, precise', parameters[:precise_discretization])
            infos.add_parameter('MP', 'minimal similarity for the 2nd pass in \'precise\' mode', parameters[:minimal_similarity])
          else
            infos.add_parameter('V', 'discretization value', parameters[:rough_discretization])
          end
          infos.background_parameter('BQ', 'background for query matrix', parameters[:query_background])
          infos.background_parameter('BC', 'background for collection', parameters[:collection_background])

          infos.add_table_parameter_without_description('motif', :name)
          infos.add_table_parameter_without_description('similarity', :similarity)
          infos.add_table_parameter_without_description('shift', :shift)
          infos.add_table_parameter_without_description('overlap', :overlap)
          infos.add_table_parameter_without_description('orientation', :orientation)
          if parameters[:precision_mode] == :precise
            infos.add_table_parameter_without_description('precise mode', :precision_mode){|precision| precision == :precise ? '*' : '.' }
          end
        }.result
      end

############################################

      def self.find_pvalue_info_string(data, parameters)
        OutputInformation.new(data) {|infos|
          infos.add_parameter('V', 'discretization value', parameters[:discretization])
          infos.background_parameter('B', 'background', parameters[:background])

          infos.add_table_parameter('T', 'threshold', :threshold)
          infos.add_table_parameter('W', 'number of recognized words', :number_of_recognized_words)  if parameters[:background] == [1,1,1,1]
          infos.add_table_parameter('P', 'P-value', :pvalue)
        }.result
      end

    end
  end
end