def read_numbers_from_file(file_path):
    """Legge una sequenza di numeri da un file e restituisce una lista di numeri interi."""
    with open(file_path, 'r') as file:
        content = file.read().strip()
    numbers = list(map(int, content.split()))
    return numbers

def compare_files(file1, file2):
    """Confronta i contenuti di due file per verificare se contengono la stessa sequenza di numeri."""
    numbers1 = read_numbers_from_file(file1)
    numbers2 = read_numbers_from_file(file2)
    
    if numbers1 == numbers2:
        print("Le sequenze di numeri nei due file sono uguali.")
    else:
        print("Le sequenze di numeri nei due file sono diverse.")

# Esempio di utilizzo
file1 = 'file1.txt'
file2 = 'file2.txt'
compare_files(file1, file2)
