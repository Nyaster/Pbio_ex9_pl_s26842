# dna_generator.py

# Cel programu:
# Program generuje losową sekwencję DNA, wstawia imię programisty w losowym miejscu tej sekwencji, zapisuje dane do pliku w formacie FASTA
# oraz oblicza statystyki występowania nukleotydów pomijając fragment z imieniem programisty.
# Może być używany do nauki analizy sekwencji genetycznych i formatowania danych biologicznych.

# KONTEKST ZASTOSOWANIA:
# Ten program może służyć jako narzędzie edukacyjne w kursach bioinformatyki, programowania lub biologii molekularnej.

import random

# Konfiguracja stałych
PROGRAMMER_NAME = "Dmytro"  # Zmień na swoje imię
FASTA_LINE_LENGTH = 60


def generate_dna_sequence(length):
    """Generuje losową sekwencję DNA o podanej długości"""
    return ''.join(random.choices(['A', 'C', 'G', 'T'], k=length))


def insert_programmer_name(sequence):
    """Wstawia imię programisty w losowym miejscu sekwencji"""
    name = PROGRAMMER_NAME
    name_len = len(name)
    seq_len = len(sequence)

    if name_len > seq_len:
        raise ValueError("Imię jest dłuższe niż sekwencja!")

    insert_pos = random.randint(0, seq_len - name_len)
    return sequence[:insert_pos] + name + sequence[insert_pos + name_len:]


def calculate_statistics(sequence, name_start, name_end):
    """Oblicza statystyki pomijając wstawione imię"""
    stats = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    total = 0

    for i, nucleotide in enumerate(sequence):
        if name_start <= i < name_end:
            continue
        if nucleotide in stats:
            stats[nucleotide] += 1
            total += 1

    return stats, total


def format_sequence(sequence):
    """Formatuje sekwencję do linii o ustalonej długości"""
    # ORIGINAL:
    # for i in range(0, len(modified_sequence), FASTA_LINE_LENGTH):
    #     fasta_file.write(modified_sequence[i:i + FASTA_LINE_LENGTH] + '\n')
    # MODIFIED (utworzono funkcję dla wielokrotnego użytku i poprawy czytelności):
    return '\n'.join(sequence[i:i + FASTA_LINE_LENGTH] for i in range(0, len(sequence), FASTA_LINE_LENGTH))


def save_to_fasta(sequence_id, description, sequence):
    """Zapisuje sekwencję do pliku w formacie FASTA"""
    with open(f"{sequence_id}.fasta", 'w') as fasta_file:
        fasta_file.write(f">{sequence_id} {description}\n")
        fasta_file.write(format_sequence(sequence) + '\n')


def display_statistics(stats, total):
    """Wyświetla statystyki nukleotydów i stosunek C+G do A+T"""
    percentages = {k: (v / total) * 100 if total > 0 else 0 for k, v in stats.items()}
    cg_ratio = (stats['C'] + stats['G']) / (stats['A'] + stats['T']) if (stats['A'] + stats['T']) > 0 else 0

    print("\nStatystyki sekwencji:")
    for nucleotide, percent in percentages.items():
        print(f"{nucleotide}: {percent:.2f}%")
    print(f"\nStosunek C+G do A+T: {cg_ratio:.2f}")


def main():
    try:
        # Pobranie danych od użytkownika
        length = int(input("Podaj długość sekwencji: "))
        if length < len(PROGRAMMER_NAME):
            # ORIGINAL:
            # sequence = generate_dna_sequence(length)
            # MODIFIED (dodano warunek sprawdzający minimalną długość sekwencji):
            raise ValueError("Długość sekwencji musi być większa niż długość imienia programisty.")

        sequence_id = input("Podaj ID sekwencji: ")
        description = input("Podaj opis sekwencji: ")

        dna_sequence = generate_dna_sequence(length)
        modified_sequence = insert_programmer_name(dna_sequence)

        name_start = modified_sequence.find(PROGRAMMER_NAME)
        name_end = name_start + len(PROGRAMMER_NAME)

        stats, total = calculate_statistics(modified_sequence, name_start, name_end)

        save_to_fasta(sequence_id, description, modified_sequence)
        display_statistics(stats, total)

    except ValueError as e:
        print(f"Błąd: {e}")
    except Exception as e:
        print(f"Nieoczekiwany błąd: {e}")


if __name__ == "__main__":
    main()
