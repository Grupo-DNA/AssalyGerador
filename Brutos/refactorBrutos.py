def refatorar_txt(input_file, output_file):
    with open(input_file, 'r') as f:
        lines = f.readlines()

    refatorado = []
    for line in lines:
        parts = line.split('\t')
        novo_line = parts[0] + '\t' + parts[1] + '\t' + parts[2] + '\t' + parts[3] + parts[4][0] + '\n'
        refatorado.append(novo_line)

    with open(output_file, 'w') as f:
        for line in refatorado:
            f.write(line)

if __name__ == "__main__":
    input_file = "CLUB003420.txt"
    output_file = "CLUB003420.txt"
    refatorar_txt(input_file, output_file)
