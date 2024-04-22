# Função para processar cada linha do arquivo
def process_line(line):
    # Separa os elementos da linha por tabulação
    elements = line.strip().split('\t')
    
    # Remove letras e '_' da terceira coluna
    elements[2] = ''.join(char for char in elements[2] if not char.isalpha() and char != '_')
    
    # Junta os elementos de volta com tabulação
    return '\t'.join(elements)

# Função principal para processar o arquivo
def process_file(input_file, output_file):
    with open(input_file, 'r') as f_input, open(output_file, 'w') as f_output:
        for line in f_input:
            processed_line = process_line(line)
            f_output.write(processed_line + '\n')

# Nome do arquivo de entrada e saída
input_file = 'CLUB003436.txt'
output_file = 'CLUB003436.txt'

# Chama a função para processar o arquivo
process_file(input_file, output_file)

print("Concluído! O arquivo foi processado com sucesso.")
