#!/bin/bash

# Script para testar todas as instâncias TSP com menos de 600 nós
# e salvar os resultados em um arquivo CSV

# Diretório das instâncias
INSTANCES_DIR="data/instances"

# Executável do programa
EXECUTABLE="build/src/project"

# Arquivo de saída
OUTPUT_FILE="results.csv"

# Verificar se o executável existe
if [ ! -f "$EXECUTABLE" ]; then
    echo "Erro: Executável não encontrado em $EXECUTABLE"
    echo "Por favor, compile o projeto primeiro com: meson compile -C build"
    exit 1
fi

# Criar cabeçalho do arquivo de resultados
echo "instance,avg_cost,avg_time_seconds" > "$OUTPUT_FILE"

# Função para extrair o número de nós de um arquivo TSP
get_dimension() {
    local file=$1
    # Procura pela linha DIMENSION no arquivo TSP
    dimension=$(grep -i "^DIMENSION" "$file" | head -n 1 | sed 's/[^0-9]*//g')
    echo "$dimension"
}

# Contador de instâncias processadas
count=0
total=0

# Contar total de instâncias < 600 nós
echo "Verificando instâncias..."
for instance in "$INSTANCES_DIR"/*.tsp; do
    [ -f "$instance" ] || continue
    dimension=$(get_dimension "$instance")
    if [ -n "$dimension" ] && [ "$dimension" -lt 600 ]; then
        total=$((total + 1))
    fi
done

echo "Encontradas $total instâncias com menos de 600 nós"
echo ""

# Processar cada instância
for instance in "$INSTANCES_DIR"/*.tsp; do
    # Verificar se o arquivo existe
    [ -f "$instance" ] || continue
    
    # Obter o nome da instância
    instance_name=$(basename "$instance" .tsp)
    
    # Obter dimensão da instância
    dimension=$(get_dimension "$instance")
    
    # Verificar se conseguiu extrair a dimensão
    if [ -z "$dimension" ]; then
        echo "Aviso: Não foi possível determinar a dimensão de $instance_name, pulando..."
        continue
    fi
    
    # Processar apenas instâncias com menos de 600 nós
    if [ "$dimension" -lt 600 ]; then
        count=$((count + 1))
        echo "[$count/$total] Processando $instance_name ($dimension nós)..."
        
        # Executar o programa e capturar a saída
        # A saída do programa é no formato: avg_cost,avg_time_seconds
        output=$("$EXECUTABLE" "$instance" 0 600)
        
        # Adicionar o nome da instância na primeira coluna
        echo "$instance_name,$output" >> "$OUTPUT_FILE"
        
        echo "  -> Concluído"
    fi
done

echo ""
echo "Processamento completo!"
echo "Resultados salvos em: $OUTPUT_FILE"
echo "Total de instâncias processadas: $count"

