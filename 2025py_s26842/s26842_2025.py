# Python script to generate nucleotide sequences in FASTA format.

# Cel programu:
# Program ma na celu generowanie losowych sekwencji DNA (składających się z nukleotydów A, C, G, T)
# o długości zdefiniowanej przez użytkownika. Umożliwia również nadanie sekwencji ID oraz opisu.
# Wynik jest zapisywany do pliku w formacie FASTA. Dodatkowo, program oblicza
# i wyświetla podstawowe statystyki wygenerowanej sekwencji, takie jak procentowa
# zawartość poszczególnych nukleotydów oraz stosunek GC. Od wersji z ulepszeniem 4,
# program potrafi również narysować wykres słupkowy tych statystyk.
# W losowym miejscu sekwencji (w wersji zapisywanej do pliku) wstawiane jest podane przez użytkownika imię,
# które nie wpływa na obliczane statystyki DNA.

# Kontekst jego zastosowania:
# Narzędzie to może być przydatne dla bioinformatyków, biologów molekularnych, studentów kierunków
# przyrodniczych lub programistów uczących się Pythona do:
# - Generowania danych testowych dla innych narzędzi bioinformatycznych.
# - Tworzenia przykładowych sekwencji do celów edukacyjnych.
# - Symulowania prostych sekwencji biologicznych.
# - Wizualizacji składu nukleotydowego.
# - Ćwiczenia w programowaniu i przetwarzaniu danych tekstowych.

import random  # Import modułu random do generowania liczb losowych i losowych wyborów
import re  # Import modułu re do pracy z wyrażeniami regularnymi (dla ulepszenia 1)

# ULEPSZENIE 4: Dodanie wizualizacji statystyk za pomocą matplotlib
# ORIGINAL:
# (brak tej funkcjonalności w oryginalnym kodzie LLM)
# MODIFIED (Dodanie opcjonalnej wizualizacji statystyk przy użyciu biblioteki matplotlib):
# Wizualizacja danych często ułatwia ich zrozumienie. Wykres słupkowy jest dobrym sposobem
# na przedstawienie procentowego udziału poszczególnych nukleotydów.
# Dodano obsługę braku biblioteki matplotlib.
MATPLOTLIB_AVAILABLE = False  # Flaga informująca, czy biblioteka matplotlib jest dostępna
try:
    import matplotlib.pyplot as plt  # Próba importu biblioteki matplotlib

    MATPLOTLIB_AVAILABLE = True  # Ustawienie flagi na True, jeśli import się powiódł
except ImportError:  # Obsługa błędu, jeśli biblioteka nie jest zainstalowana
    print("Biblioteka matplotlib nie jest zainstalowana. Wykresy nie będą generowane.")  # Informacja dla użytkownika
    print("Aby zainstalować, użyj: pip install matplotlib")


# --- Funkcje oryginalne i zmodyfikowane ---

def generuj_sekwencje_dna(dlugosc: int) -> str:  # Definicja funkcji generującej sekwencję DNA
    """Generuje losową sekwencję DNA o podanej długości."""  # Docstring funkcji
    if dlugosc < 0:  # Sprawdzenie, czy podana długość nie jest ujemna
        raise ValueError("Długość sekwencji nie może być ujemna.")  # Podniesienie błędu, jeśli warunek jest spełniony
    # ORIGINAL:
    # return "".join(random.choice("ACGT") for _ in range(dlugosc))
    # MODIFIED (Ulepszenie poza wymaganymi 3: dodanie obsługi długości 0 w bardziej явny sposób, choć poprzednia wersja też by zadziałała poprawnie zwracając pusty string):
    # Chociaż oryginalny kod działałby poprawnie dla długości 0, zwracając pusty ciąg,
    # ta modyfikacja czyni obsługę tego przypadku bardziej jawną i czytelną.
    if dlugosc == 0:  # Jeśli długość to 0
        return ""  # Zwróć pusty string
    return "".join(random.choice("ACGT") for _ in range(
        dlugosc))  # Generowanie sekwencji: pętla przez zadaną długość, w każdej iteracji losowy wybór z "ACGT" i połączenie w string


def oblicz_statystyki(sekwencja: str) -> dict:  # Definicja funkcji obliczającej statystyki sekwencji
    """Oblicza statystyki dla podanej sekwencji DNA, w tym częstości dinukleotydów."""  # Docstring funkcji
    # Inicjalizacja słownika na statystyki mononukleotydów
    statystyki_mono_counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}  # Zmieniono nazwę dla jasności (surowe liczby)
    dlugosc_sek = len(sekwencja)  # Pobranie długości analizowanej sekwencji

    if dlugosc_sek == 0:  # Jeśli sekwencja jest pusta
        return {  # Zwróć słownik z zerowymi wartościami dla statystyk
            'A_proc': 0.0, 'C_proc': 0.0, 'G_proc': 0.0, 'T_proc': 0.0,
            'A_count': 0, 'C_count': 0, 'G_count': 0, 'T_count': 0,  # Dodano surowe liczby
            'CG_stosunek': 0.0,
            'dinukleotydy': {}  # Pusty słownik dla dinukleotydów
        }

    for nukleotyd in sekwencja:  # Pętla przez każdy nukleotyd w sekwencji
        if nukleotyd in statystyki_mono_counts:  # Sprawdzenie, czy nukleotyd jest jednym z A, C, G, T
            statystyki_mono_counts[nukleotyd] += 1  # Inkrementacja licznika dla danego nukleotydu

    # Obliczanie procentowej zawartości mononukleotydów
    statystyki_proc = {
        'A_proc': (statystyki_mono_counts['A'] / dlugosc_sek) * 100,  # Procent A
        'C_proc': (statystyki_mono_counts['C'] / dlugosc_sek) * 100,  # Procent C
        'G_proc': (statystyki_mono_counts['G'] / dlugosc_sek) * 100,  # Procent G
        'T_proc': (statystyki_mono_counts['T'] / dlugosc_sek) * 100,  # Procent T
        'A_count': statystyki_mono_counts['A'],  # Surowa liczba A
        'C_count': statystyki_mono_counts['C'],  # Surowa liczba C
        'G_count': statystyki_mono_counts['G'],  # Surowa liczba G
        'T_count': statystyki_mono_counts['T'],  # Surowa liczba T
    }

    liczba_cg = statystyki_mono_counts['C'] + statystyki_mono_counts['G']  # Suma liczby C i G
    statystyki_proc['CG_stosunek'] = (liczba_cg / dlugosc_sek) * 100  # Procentowa zawartość GC

    # ULEPSZENIE 3: Dodanie obliczania częstości dinukleotydów
    # ORIGINAL:
    # (brak tej funkcjonalności w oryginalnym kodzie LLM, statystyki_proc zawierały tylko A,C,G,T proc i CG_stosunek)
    # MODIFIED (Dodanie obliczania częstości dinukleotydów dla pełniejszej analizy sekwencji):
    # Biologicznie, częstości dinukleotydów (np. CpG) są ważne. Ta zmiana rozszerza analizę.
    statystyki_dinukleotydow = {}  # Inicjalizacja słownika na częstości dinukleotydów
    if dlugosc_sek > 1:  # Obliczanie dinukleotydów ma sens tylko jeśli sekwencja ma co najmniej 2 nukleotydy
        for i in range(dlugosc_sek - 1):  # Pętla od pierwszego do przedostatniego nukleotydu
            dinukleotyd = sekwencja[i:i + 2]  # Pobranie pary nukleotydów (dinukleotydu)
            statystyki_dinukleotydow[dinukleotyd] = statystyki_dinukleotydow.get(dinukleotyd,
                                                                                 0) + 1  # Inkrementacja licznika dla danego dinukleotydu
    statystyki_proc['dinukleotydy'] = statystyki_dinukleotydow  # Dodanie statystyk dinukleotydów do wynikowego słownika

    return statystyki_proc  # Zwrócenie słownika ze wszystkimi obliczonymi statystykami


def wstaw_imie(sekwencja: str, imie: str) -> str:  # Definicja funkcji wstawiającej imię do sekwencji
    """Wstawia imię w losowe miejsce sekwencji."""  # Docstring funkcji
    if not sekwencja:  # Jeśli sekwencja jest pusta
        return imie  # Zwróć samo imię (lub pusty string jeśli imię też puste)
    if not imie:  # Jeśli imię jest puste
        return sekwencja  # Zwróć oryginalną, nienaruszoną sekwencję

    pozycja = random.randint(0,
                             len(sekwencja))  # Wybór losowej pozycji do wstawienia imienia (może być na początku, w środku lub na końcu)
    return sekwencja[:pozycja] + imie + sekwencja[
                                        pozycja:]  # Zwrócenie nowej sekwencji: część przed pozycją + imię + część po pozycji


# ULEPSZENIE 2: Modyfikacja funkcji zapisującej do FASTA w celu zawijania długich linii sekwencji
def zapisz_do_fasta(nazwa_pliku: str, id_sekwencji: str, opis: str, sekwencja_z_imieniem: str,
                    szerokosc_linii: int = 70):  # Definicja funkcji zapisującej do pliku FASTA
    """Zapisuje sekwencję do pliku w formacie FASTA, z opcjonalnym zawijaniem linii."""  # Docstring funkcji
    with open(nazwa_pliku, 'w') as plik:  # Otwarcie pliku w trybie zapisu (tryb 'w' nadpisuje plik, jeśli istnieje)
        plik.write(f">{id_sekwencji} {opis}\n")  # Zapis nagłówka FASTA (linia zaczynająca się od '>')
        # ORIGINAL:
        # plik.write(f"{sekwencja_z_imieniem}\n")
        # MODIFIED (Dodanie zawijania linii sekwencji dla lepszej czytelności plików FASTA zgodnie ze standardem):
        # Standard FASTA często zakłada zawijanie linii sekwencji (np. co 60-80 znaków), co ułatwia przeglądanie.
        if szerokosc_linii <= 0:  # Jeśli podano niepoprawną szerokość linii (lub 0 dla braku zawijania)
            plik.write(f"{sekwencja_z_imieniem}\n")  # Zapisz całą sekwencję w jednej linii
        else:
            for i in range(0, len(sekwencja_z_imieniem),
                           szerokosc_linii):  # Pętla przez sekwencję z krokiem równym szerokości linii
                plik.write(sekwencja_z_imieniem[
                           i:i + szerokosc_linii] + "\n")  # Zapis fragmentu sekwencji i przejście do nowej linii


# ULEPSZENIE 1: Funkcja do oczyszczania ID sekwencji na potrzeby nazwy pliku
def oczysc_id_dla_nazwy_pliku(id_sekwencji: str) -> str:  # Definicja funkcji oczyszczającej ID
    """Oczyszcza ID sekwencji z niebezpiecznych znaków, aby mogło być użyte jako nazwa pliku."""  # Docstring funkcji
    # ORIGINAL:
    # W oryginalnym kodzie LLM było tylko: nazwa_pliku_fasta = f"{id_sekwencji.replace(' ', '_')}.fasta"
    # To obsługiwało tylko spacje.
    # MODIFIED (Zastąpienie niebezpiecznych znaków w ID sekwencji, aby uniknąć błędów przy tworzeniu pliku):
    # Nazwy plików nie mogą zawierać pewnych znaków (np. / \ : * ? " < > |). Ta funkcja usuwa je lub zastępuje.
    id_oczyszczone = re.sub(r'[\\/*?:"<>|\s]', '_',
                            id_sekwencji)  # Użycie wyrażenia regularnego do zastąpienia spacji i innych niebezpiecznych znaków podkreślnikiem
    id_oczyszczone = re.sub(r'_+', '_', id_oczyszczone)  # Zastąpienie wielokrotnych podkreślników pojedynczym
    if not id_oczyszczone:  # Jeśli po oczyszczeniu ID jest puste
        return "default_seq_id"  # Zwróć domyślną nazwę
    return id_oczyszczone  # Zwróć oczyszczone ID


# ULEPSZENIE 4 (ciąg dalszy): Funkcja rysująca wykres statystyk
def rysuj_statystyki_wykres(statystyki: dict, id_sekwencji: str,
                            dlugosc_oryginalna: int):  # Definicja funkcji rysującej wykres
    """Rysuje wykres słupkowy procentowej zawartości nukleotydów."""  # Docstring funkcji
    if not MATPLOTLIB_AVAILABLE:  # Sprawdzenie, czy biblioteka matplotlib jest dostępna
        print(
            "Wykres nie może być wygenerowany, ponieważ biblioteka matplotlib nie jest dostępna.")  # Informacja dla użytkownika
        return  # Zakończenie funkcji, jeśli biblioteka nie jest dostępna

    labels = ['A', 'C', 'G', 'T']  # Etykiety dla osi X (nukleotydy)
    procenty = [statystyki['A_proc'], statystyki['C_proc'], statystyki['G_proc'],
                statystyki['T_proc']]  # Wartości procentowe dla każdego nukleotydu

    x = range(len(labels))  # Pozycje dla słupków na osi X

    fig, ax = plt.subplots()  # Utworzenie figury i osi wykresu
    rects = ax.bar(x, procenty, label='Procenty')  # Narysowanie słupków

    # Dodanie etykiet, tytułu i legendy
    ax.set_ylabel('Procentowa zawartość (%)')  # Etykieta osi Y
    ax.set_xlabel('Nukleotyd')  # Etykieta osi X
    ax.set_title(
        f'Statystyki nukleotydów dla sekwencji: {id_sekwencji}\nOryginalna długość: {dlugosc_oryginalna} bp, %GC: {statystyki["CG_stosunek"]:.1f}%')  # Tytuł wykresu
    ax.set_xticks(x)  # Ustawienie pozycji etykiet na osi X
    ax.set_xticklabels(labels)  # Ustawienie etykiet na osi X
    # ax.legend() # Legenda (opcjonalna, przy jednym zestawie danych może być zbędna)

    # Dodanie wartości procentowych nad słupkami
    for rect in rects:  # Pętla przez każdy słupek
        height = rect.get_height()  # Pobranie wysokości słupka
        ax.annotate(f'{height:.1f}%',  # Tekst adnotacji (wartość procentowa)
                    xy=(rect.get_x() + rect.get_width() / 2, height),  # Pozycja XY adnotacji (środek góry słupka)
                    xytext=(0, 3),  # Przesunięcie tekstu w pionie (3 punkty nad słupkiem)
                    textcoords="offset points",  # System współrzędnych dla przesunięcia
                    ha='center', va='bottom')  # Wyrównanie tekstu

    fig.tight_layout()  # Dopasowanie układu, aby elementy się nie nakładały
    plt.show()  # Wyświetlenie wykresu


def main():  # Główna funkcja programu, sterująca jego wykonaniem
    """Główna funkcja programu."""  # Docstring funkcji
    print("Generator sekwencji DNA w formacie FASTA")  # Wyświetlenie tytułu programu
    print("-----------------------------------------")  # Wyświetlenie separatora

    while True:  # Pętla nieskończona do pobierania długości sekwencji (przerywana przez break)
        try:  # Blok try-except do obsługi błędów przy konwersji inputu na liczbę
            dlugosc = int(
                input("Podaj długość sekwencji DNA (np. 100): "))  # Pobranie długości sekwencji od użytkownika
            if dlugosc < 0:  # Sprawdzenie, czy długość nie jest ujemna
                print("Długość nie może być ujemna. Spróbuj ponownie.")  # Komunikat o błędzie
                continue  # Powrót na początek pętli while
            break  # Wyjście z pętli, jeśli długość jest poprawna
        except ValueError:  # Obsługa błędu, jeśli użytkownik podał wartość, której nie da się skonwertować na int
            print("Nieprawidłowa wartość. Długość musi być liczbą całkowitą.")  # Komunikat o błędzie

    id_sekwencji_uzytkownika = ""  # Inicjalizacja zmiennej na ID sekwencji
    while not id_sekwencji_uzytkownika.strip():  # Pętla dopóki ID podane przez użytkownika nie będzie zawierało znaków innych niż białe
        id_sekwencji_uzytkownika = input("Podaj ID sekwencji (np. Seq1): ")  # Pobranie ID sekwencji od użytkownika
        if not id_sekwencji_uzytkownika.strip():  # Sprawdzenie, czy ID po usunięciu białych znaków z początku/końca jest puste
            print("ID sekwencji nie może być puste. Spróbuj ponownie.")  # Komunikat o błędzie

    opis_sekwencji = input("Podaj opis sekwencji (np. Losowa sekwencja testowa): ")  # Pobranie opisu sekwencji
    imie_ai = input("Podaj swoje imię (np. Gemini) do wstawienia w sekwencję: ")  # Pobranie imienia do wstawienia

    # 1. Generuj czystą sekwencję DNA
    czysta_sekwencja_dna = generuj_sekwencje_dna(dlugosc)  # Wywołanie funkcji generującej sekwencję DNA

    # 2. Oblicz statystyki na podstawie czystej sekwencji DNA
    staty = oblicz_statystyki(czysta_sekwencja_dna)  # Wywołanie funkcji obliczającej statystyki

    # 3. Wstaw imię do sekwencji (ta wersja będzie zapisana do pliku)
    sekwencja_do_zapisu = wstaw_imie(czysta_sekwencja_dna, imie_ai)  # Wywołanie funkcji wstawiającej imię

    # 4. Przygotuj nazwę pliku i zapisz do pliku
    # ORIGINAL (w kontekście tworzenia nazwy pliku):
    # nazwa_pliku_fasta = f"{id_sekwencji.replace(' ', '_')}.fasta"
    # MODIFIED (Użycie funkcji oczysc_id_dla_nazwy_pliku dla bezpieczniejszej nazwy pliku - implementacja ULEPSZENIA 1):
    # Zastosowanie funkcji oczyszczającej ID zapewnia, że nazwa pliku będzie poprawna w systemie plików.
    oczyszczone_id = oczysc_id_dla_nazwy_pliku(id_sekwencji_uzytkownika)  # Oczyszczenie ID na potrzeby nazwy pliku
    nazwa_pliku_fasta = f"{oczyszczone_id}.fasta"  # Utworzenie pełnej nazwy pliku z rozszerzeniem .fasta

    zapisz_do_fasta(nazwa_pliku_fasta, id_sekwencji_uzytkownika, opis_sekwencji,
                    sekwencja_do_zapisu)  # Wywołanie funkcji zapisującej do pliku FASTA
    print(f"\nSekwencja została zapisana do pliku {nazwa_pliku_fasta}")  # Komunikat o sukcesie zapisu

    # 5. Wyświetl statystyki (obliczone na podstawie czystej sekwencji)
    print("Statystyki sekwencji (na podstawie oryginalnej sekwencji DNA):")  # Nagłówek dla statystyk
    if dlugosc > 0:  # Jeśli długość sekwencji jest większa od zera
        print(f"  Długość oryginalnej sekwencji DNA: {dlugosc}")  # Wyświetlenie długości
        print(f"  A: {staty['A_proc']:.1f}% ({staty['A_count']})")  # Wyświetlenie procentu i liczby A
        print(f"  C: {staty['C_proc']:.1f}% ({staty['C_count']})")  # Wyświetlenie procentu i liczby C
        print(f"  G: {staty['G_proc']:.1f}% ({staty['G_count']})")  # Wyświetlenie procentu i liczby G
        print(f"  T: {staty['T_proc']:.1f}% ({staty['T_count']})")  # Wyświetlenie procentu i liczby T
        print(f"  %CG: {staty['CG_stosunek']:.1f}%")  # Wyświetlenie procentu GC
        if staty['dinukleotydy']:  # Jeśli są dostępne statystyki dinukleotydów
            print("  Częstości dinukleotydów:")  # Nagłówek dla dinukleotydów
            # Sortowanie dinukleotydów alfabetycznie dla spójnego wyświetlania
            for dinukl, liczba in sorted(staty['dinukleotydy'].items()):  # Pętla przez posortowane dinukleotydy
                print(f"    {dinukl}: {liczba}")  # Wyświetlenie dinukleotydu i jego liczby
        else:
            print("  (Brak dinukleotydów do wyświetlenia - sekwencja zbyt krótka)")  # Komunikat o braku dinukleotydów
    else:  # Jeśli długość sekwencji wynosi 0
        print("  Długość oryginalnej sekwencji DNA: 0")  # Informacja o zerowej długości
        print("  A: 0.0% (0)")  # Statystyki dla zerowej długości
        print("  C: 0.0% (0)")
        print("  G: 0.0% (0)")
        print("  T: 0.0% (0)")
        print("  %CG: 0.0%")
        print("  (Brak dinukleotydów do wyświetlenia - sekwencja o długości 0)")

    # 6. Narysuj wykres statystyk (ULEPSZENIE 4)
    if dlugosc > 0 and MATPLOTLIB_AVAILABLE:  # Rysuj wykres tylko jeśli są dane i biblioteka jest dostępna
        print("\nGenerowanie wykresu statystyk...")  # Informacja dla użytkownika
        rysuj_statystyki_wykres(staty, id_sekwencji_uzytkownika, dlugosc)  # Wywołanie funkcji rysującej wykres
    elif dlugosc == 0 and MATPLOTLIB_AVAILABLE:
        print("\nNie można wygenerować wykresu dla sekwencji o długości 0.")


if __name__ == "__main__":  # Standardowy idiom w Pythonie: kod w tym bloku wykona się tylko, gdy plik jest uruchamiany bezpośrednio (nie importowany jako moduł)
    main()  # Wywołanie głównej funkcji programu