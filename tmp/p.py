import numpy as np

# Parametri della distribuzione
media = 207.756210  
varianza =  2.2
deviazione_standard = np.sqrt(varianza)
n_tempi = 9

# Generazione dei tempi di esecuzione
tempi_esecuzione = np.random.normal(media, deviazione_standard, n_tempi)
tempi_esecuzione = np.round(tempi_esecuzione, 6)
print(tempi_esecuzione)
