# Função para processar os dados dos SNPs e armazenar os valores dos SNPs de interesse
def processar_snps(snps_data):
    apoe_snps = {'rs429358': None, 'rs7412': None}
    for snp in snps_data:
        parts = snp.split('\t')
        snp_id = parts[0]
        snp_value = parts[1]
        if snp_id in apoe_snps:
            apoe_snps[snp_id] = snp_value
    return apoe_snps

# Função para determinar a variante APOE com base nos valores dos SNPs
def determinar_apoe_variantes(rs429358, rs7412):
    if rs429358 == 'T' and rs7412 == 'T':
        return 'ε2/ε2'
    elif rs429358 == 'T' and rs7412 == 'C':
        return 'ε2/ε3'
    elif rs429358 == 'C' and rs7412 == 'T':
        return 'ε2/ε4'
    elif rs429358 == 'C' and rs7412 == 'C':
        return 'ε4/ε4'
    elif rs429358 == 'T' and rs7412 == 'C':
        return 'ε3/ε3'
    return 'Desconhecido'

# Função principal que coordena o processamento e a determinação da variante APOE
def determinar_variante_apoe(snps_data):
    apoe_snps = processar_snps(snps_data)
    
    if apoe_snps['rs429358'] and apoe_snps['rs7412']:
        rs429358 = apoe_snps['rs429358'][0]
        rs7412 = apoe_snps['rs7412'][0]
        apoe_variantes = determinar_apoe_variantes(rs429358, rs7412)
        return f'Variante APOE do paciente: {apoe_variantes}'
    else:
        return 'Dados dos SNPs rs429358 e rs7412 não encontrados nos dados do paciente'