"""
Inorganic Chemistry API Router — расширенный функционал для неорганической химии
"""

from fastapi import APIRouter, HTTPException, Query, Depends
from pydantic import BaseModel, Field
from typing import List, Optional, Dict, Any
import urllib.parse
import urllib.request
import json
import math

from database import get_db
from Molecule import Molecule
from i18n import get_text, get_lang

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors

router = APIRouter(prefix="/inorg/api", tags=["Inorganic Chemistry API"])

PROXY_ADDRESS = "127.0.0.1:9050"


def set_proxy(proxy: str):
    global PROXY_ADDRESS
    PROXY_ADDRESS = proxy


# ─────────────────────────────────────────────────────────────────────────────
# Константы
# ─────────────────────────────────────────────────────────────────────────────

# Электроотрицательность Полинга
ELECTRONEGATIVITY = {
    "H": 2.20, "He": 0.00,
    "Li": 0.98, "Be": 1.57, "B": 2.04, "C": 2.55, "N": 3.04, "O": 3.44, "F": 3.98, "Ne": 0.00,
    "Na": 0.93, "Mg": 1.31, "Al": 1.61, "Si": 1.90, "P": 2.19, "S": 2.58, "Cl": 3.16, "Ar": 0.00,
    "K": 0.82, "Ca": 1.00, "Sc": 1.36, "Ti": 1.54, "V": 1.63, "Cr": 1.66, "Mn": 1.55, "Fe": 1.83,
    "Co": 1.88, "Ni": 1.91, "Cu": 1.90, "Zn": 1.65, "Ga": 1.81, "Ge": 2.01, "As": 2.18, "Se": 2.55,
    "Br": 2.96, "Kr": 3.00,
    "Rb": 0.82, "Sr": 0.95, "Y": 1.22, "Zr": 1.33, "Nb": 1.60, "Mo": 2.16, "Tc": 1.90, "Ru": 2.20,
    "Rh": 2.28, "Pd": 2.20, "Ag": 1.93, "Cd": 1.69, "In": 1.78, "Sn": 1.96, "Sb": 2.05, "Te": 2.10,
    "I": 2.66, "Xe": 2.60,
    "Cs": 0.79, "Ba": 0.89, "La": 1.10, "Hf": 1.30, "Ta": 1.50, "W": 2.36, "Re": 1.90, "Os": 2.20,
    "Ir": 2.20, "Pt": 2.28, "Au": 2.54, "Hg": 2.00, "Tl": 1.62, "Pb": 2.33, "Bi": 2.02, "Po": 2.00,
    "At": 2.20, "Rn": 0.00,
    "Fr": 0.70, "Ra": 0.90, "Ac": 1.10, "Th": 1.30, "Pa": 1.50, "U": 1.38, "Np": 1.36, "Pu": 1.28,
    "Am": 1.30, "Cm": 1.30, "Bk": 1.30, "Cf": 1.30, "Es": 1.30, "Fm": 1.30, "Md": 1.30, "No": 1.30,
}

# Атомные массы
ATOMIC_MASSES = {
    "H": 1.008, "He": 4.003, "Li": 6.941, "Be": 9.012, "B": 10.81, "C": 12.011, "N": 14.007,
    "O": 15.999, "F": 18.998, "Ne": 20.180, "Na": 22.990, "Mg": 24.305, "Al": 26.982,
    "Si": 28.086, "P": 30.974, "S": 32.065, "Cl": 35.453, "Ar": 39.948, "K": 39.098,
    "Ca": 40.078, "Sc": 44.956, "Ti": 47.867, "V": 50.942, "Cr": 51.996, "Mn": 54.938,
    "Fe": 55.845, "Co": 58.933, "Ni": 58.693, "Cu": 63.546, "Zn": 65.380, "Ga": 69.723,
    "Ge": 72.630, "As": 74.922, "Se": 78.971, "Br": 79.904, "Kr": 83.798, "Rb": 85.468,
    "Sr": 87.620, "Y": 88.906, "Zr": 91.224, "Nb": 92.906, "Mo": 95.950, "Tc": 98.000,
    "Ru": 101.070, "Rh": 102.906, "Pd": 106.420, "Ag": 107.868, "Cd": 112.414, "In": 114.818,
    "Sn": 118.710, "Sb": 121.760, "Te": 127.600, "I": 126.904, "Xe": 131.293, "Cs": 132.905,
    "Ba": 137.327, "La": 138.905, "Hf": 178.490, "Ta": 180.948, "W": 183.840, "Re": 186.207,
    "Os": 190.230, "Ir": 192.217, "Pt": 195.084, "Au": 196.967, "Hg": 200.592, "Tl": 204.383,
    "Pb": 207.200, "Bi": 208.980, "Po": 209.000, "At": 210.000, "Rn": 222.000, "Fr": 223.000,
    "Ra": 226.000, "Ac": 227.000, "Th": 232.038, "Pa": 231.036, "U": 238.029,
}

# Типичные степени окисления
OXIDATION_STATES = {
    "H": [-1, 1], "Li": [1], "Be": [2], "B": [3], "C": [-4, -3, -2, -1, 0, 1, 2, 3, 4],
    "N": [-3, -2, -1, 0, 1, 2, 3, 4, 5], "O": [-2, -1, 0, 1, 2], "F": [-1],
    "Na": [1], "Mg": [2], "Al": [3], "Si": [-4, 4], "P": [-3, 3, 5], "S": [-2, 2, 4, 6],
    "Cl": [-1, 1, 3, 5, 7], "K": [1], "Ca": [2], "Sc": [3], "Ti": [2, 3, 4], "V": [2, 3, 4, 5],
    "Cr": [2, 3, 6], "Mn": [2, 3, 4, 6, 7], "Fe": [2, 3], "Co": [2, 3], "Ni": [2, 3],
    "Cu": [1, 2], "Zn": [2], "Ga": [3], "Ge": [2, 4], "As": [-3, 3, 5], "Se": [-2, 4, 6],
    "Br": [-1, 1, 3, 5, 7], "Rb": [1], "Sr": [2], "Y": [3], "Zr": [4], "Nb": [3, 5],
    "Mo": [2, 3, 4, 5, 6], "Tc": [4, 7], "Ru": [2, 3, 4, 6, 8], "Rh": [2, 3, 4],
    "Pd": [2, 4], "Ag": [1], "Cd": [2], "In": [1, 3], "Sn": [2, 4], "Sb": [-3, 3, 5],
    "Te": [-2, 4, 6], "I": [-1, 1, 3, 5, 7], "Cs": [1], "Ba": [2], "La": [3],
    "Hf": [4], "Ta": [5], "W": [2, 3, 4, 5, 6], "Re": [2, 3, 4, 6, 7], "Os": [2, 3, 4, 6, 8],
    "Ir": [2, 3, 4, 6], "Pt": [2, 4], "Au": [1, 3], "Hg": [1, 2], "Tl": [1, 3],
    "Pb": [2, 4], "Bi": [3, 5], "Po": [2, 4], "At": [-1, 1, 3, 5, 7],
    "Fr": [1], "Ra": [2], "Ac": [3], "Th": [4], "Pa": [4, 5], "U": [3, 4, 5, 6],
}

# Таблица растворимости (упрощенная)
# Ключи: катион, анион -> растворимость
SOLUBILITY_TABLE = {
    # Нитраты, ацетаты, хлораты, перхлораты - всегда растворимы
    "NO3": "soluble", "CH3COO": "soluble", "ClO3": "soluble", "ClO4": "soluble",
    # Хлориды, бромиды, иодиды
    "Cl": {"Ag": "insoluble", "Pb": "slightly", "Hg2": "insoluble", "Cu": "insoluble", "Tl": "insoluble"},
    "Br": {"Ag": "insoluble", "Pb": "slightly", "Hg2": "insoluble", "Cu": "insoluble", "Tl": "insoluble"},
    "I": {"Ag": "insoluble", "Pb": "insoluble", "Hg2": "insoluble", "Cu": "insoluble", "Tl": "insoluble", "Hg": "insoluble"},
    # Сульфаты
    "SO4": {"Ba": "insoluble", "Pb": "insoluble", "Ca": "slightly", "Sr": "insoluble", "Ag": "slightly", "Hg2": "insoluble"},
    # Сульфиды
    "S": {"NH4": "soluble", "alkali": "soluble", "alkaline_earth": "soluble", "Mg": "slightly", "Ca": "slightly", "Sr": "slightly", "Ba": "slightly"},
    # Карбонаты
    "CO3": {"NH4": "soluble", "alkali": "soluble", "Mg": "slightly", "Ca": "insoluble", "Sr": "insoluble", "Ba": "insoluble", "Fe": "insoluble", "Cu": "insoluble", "Zn": "insoluble", "Ag": "insoluble"},
    # Фосфаты
    "PO4": {"NH4": "soluble", "alkali": "soluble", "Mg": "slightly", "Ca": "insoluble", "Sr": "insoluble", "Ba": "insoluble", "Fe": "insoluble", "Al": "insoluble"},
    # Гидроксиды
    "OH": {"NH4": "soluble", "alkali": "soluble", "Ba": "soluble", "Sr": "soluble", "Ca": "slightly", "Mg": "insoluble", "Al": "insoluble", "Fe": "insoluble", "Cu": "insoluble", "Zn": "insoluble", "Ag": "insoluble"},
    # Оксалаты
    "C2O4": {"NH4": "soluble", "alkali": "soluble", "Mg": "slightly", "Ca": "insoluble", "Sr": "insoluble", "Ba": "insoluble", "Fe": "insoluble", "Cu": "insoluble", "Zn": "insoluble"},
    # Хроматы
    "CrO4": {"NH4": "soluble", "alkali": "soluble", "Mg": "soluble", "Ca": "insoluble", "Sr": "insoluble", "Ba": "insoluble", "Pb": "insoluble", "Ag": "insoluble"},
    # Перманганаты
    "MnO4": {"NH4": "soluble", "alkali": "soluble"},
}

# Щелочные металлы
ALKALI_METALS = {"Li", "Na", "K", "Rb", "Cs", "Fr"}
# Щёлочноземельные металлы
ALKALINE_EARTH = {"Be", "Mg", "Ca", "Sr", "Ba", "Ra"}


# ─────────────────────────────────────────────────────────────────────────────
# Вспомогательные функции
# ─────────────────────────────────────────────────────────────────────────────

def _get_json(url: str, timeout: int = 10) -> dict:
    """Выполняет HTTP GET запрос с поддержкой SOCKS5 прокси."""
    import socks
    from sockshandler import SocksiPyHandler

    for attempt_proxy in [True, False]:
        try:
            if attempt_proxy and PROXY_ADDRESS.strip():
                proxy_host, proxy_port = PROXY_ADDRESS.split(":")
                proxy_handler = SocksiPyHandler(socks.SOCKS5, proxy_host, int(proxy_port))
                opener = urllib.request.build_opener(proxy_handler)
            else:
                opener = urllib.request.build_opener()

            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
            with opener.open(req, timeout=timeout) as resp:
                return json.loads(resp.read().decode("utf-8"))
        except Exception:
            if attempt_proxy:
                continue
            raise

    raise HTTPException(status_code=502, detail="External API unavailable")


def _parse_formula(formula: str) -> Dict[str, int]:
    """
    Парсит химическую формулу и возвращает словарь {элемент: количество}.
    Поддерживает скобки и коэффициенты.
    """
    import re
    tokens = re.findall(r'([A-Z][a-z]?)(\d*)|(\()|(\))(\d*)', formula)
    stack = [{}]
    for elem, count, open_paren, close_paren, close_count in tokens:
        if elem:
            n = int(count) if count else 1
            stack[-1][elem] = stack[-1].get(elem, 0) + n
        elif open_paren:
            stack.append({})
        elif close_paren:
            n = int(close_count) if close_count else 1
            top = stack.pop()
            for e, c in top.items():
                stack[-1][e] = stack[-1].get(e, 0) + c * n
    return stack[0]


def _get_element_info(symbol: str) -> dict:
    """Получает информацию об элементе из PubChem PUG REST API."""
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/element/symbol/{symbol}/JSON"
        data = _get_json(url)
        if "Table" in data and "Row" in data["Table"]:
            row = data["Table"]["Row"][0]
            return {
                "symbol": symbol,
                "atomic_number": row.get("Cell", [None])[0],
                "atomic_mass": row.get("Cell", [None, None])[1],
            }
    except Exception:
        pass
    # Fallback на локальные данные
    return {
        "symbol": symbol,
        "atomic_number": None,
        "atomic_mass": ATOMIC_MASSES.get(symbol),
    }


# ─────────────────────────────────────────────────────────────────────────────
# Эндпоинты
# ─────────────────────────────────────────────────────────────────────────────

@router.get("/", include_in_schema=False)
async def inorg_index(lang: str = Depends(get_lang)):
    """Индексная страница Inorganic Chemistry API."""
    return {
        "endpoints": {
            "GET /inorg/api/parse_formula": "Парсинг химической формулы",
            "GET /inorg/api/moles": "Расчёт молей и массы",
            "GET /inorg/api/electronegativity": "Электроотрицательность и тип связи",
            "GET /inorg/api/oxidation_states": "Степени окисления",
            "GET /inorg/api/solubility": "Проверка растворимости",
            "GET /inorg/api/redox": "Окислительно-восстановительные реакции",
            "GET /inorg/api/stoichiometry": "Стехиометрия реакции",
            "GET /inorg/api/orbitals": "Электронная конфигурация и орбитали",
            "GET /inorg/api/element_info": "Информация об элементе",
        }
    }


@router.get("/parse_formula", summary="Парсинг химической формулы")
def parse_formula(
    formula: str = Query(..., description="Химическая формула, например H2SO4"),
    lang: str = Depends(get_lang),
):
    """
    Парсит химическую формулу и возвращает:
    - Список элементов с количествами
    - Молекулярную массу
    - Массовые доли элементов
    """
    try:
        elements = _parse_formula(formula)
        if not elements:
            raise ValueError("Не удалось разобрать формулу")

        # Вычисляем молекулярную массу
        total_mass = 0.0
        element_details = []
        for elem, count in elements.items():
            mass = ATOMIC_MASSES.get(elem)
            if mass is None:
                raise ValueError(f"Неизвестный элемент: {elem}")
            total_mass += mass * count
            element_details.append({
                "element": elem,
                "count": count,
                "atomic_mass": mass,
                "total_mass": round(mass * count, 4),
            })

        # Массовые доли
        for ed in element_details:
            ed["mass_percent"] = round(ed["total_mass"] / total_mass * 100, 2)

        return {
            "formula": formula,
            "elements": element_details,
            "molecular_mass": round(total_mass, 4),
            "total_atoms": sum(elements.values()),
        }
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


@router.get("/moles", summary="Расчёт молей и массы")
def calculate_moles(
    formula: str = Query(..., description="Химическая формула"),
    mass: float = Query(..., description="Масса в граммах", gt=0),
    lang: str = Depends(get_lang),
):
    """
    Рассчитывает:
    - Количество молей по массе
    - Массу по количеству молей
    - Число молекул (число Авогадро)
    """
    try:
        elements = _parse_formula(formula)
        if not elements:
            raise ValueError("Не удалось разобрать формулу")

        # Молекулярная масса
        total_mass = sum(ATOMIC_MASSES.get(elem, 0) * count for elem, count in elements.items())
        if total_mass == 0:
            raise ValueError("Не удалось вычислить молекулярную массу")

        moles = mass / total_mass
        avogadro = 6.02214076e23
        molecules = moles * avogadro

        return {
            "formula": formula,
            "mass_g": mass,
            "molecular_mass": round(total_mass, 4),
            "moles": round(moles, 6),
            "molecules": molecules,
            "atoms": sum(elements.values()),
        }
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


@router.get("/electronegativity", summary="Электроотрицательность и тип связи")
def get_electronegativity(
    element1: str = Query(..., description="Первый элемент (например, Na)"),
    element2: str = Query(..., description="Второй элемент (например, Cl)"),
    lang: str = Depends(get_lang),
):
    """
    Определяет тип связи между двумя элементами по разнице электроотрицательностей:
    - ΔEN < 0.5: неполярная ковалентная
    - 0.5 ≤ ΔEN < 1.7: полярная ковалентная
    - ΔEN ≥ 1.7: ионная
    """
    try:
        en1 = ELECTRONEGATIVITY.get(element1)
        en2 = ELECTRONEGATIVITY.get(element2)

        if en1 is None or en2 is None:
            raise ValueError(f"Неизвестный элемент: {element1 if en1 is None else element2}")

        delta = abs(en1 - en2)

        if delta < 0.5:
            bond_type = "nonpolar_covalent"
            bond_type_ru = "неполярная ковалентная"
        elif delta < 1.7:
            bond_type = "polar_covalent"
            bond_type_ru = "полярная ковалентная"
        else:
            bond_type = "ionic"
            bond_type_ru = "ионная"

        return {
            "element1": element1,
            "element2": element2,
            "electronegativity1": en1,
            "electronegativity2": en2,
            "delta": round(delta, 2),
            "bond_type": bond_type,
            "bond_type_ru": bond_type_ru,
        }
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


@router.get("/oxidation_states", summary="Степени окисления")
def get_oxidation_states(
    formula: str = Query(..., description="Химическая формула, например Fe2O3"),
    lang: str = Depends(get_lang),
):
    """
    Определяет степени окисления элементов в соединении.
    Использует правила: O = -2, H = +1, щелочные = +1, галогены = -1 (обычно).
    """
    try:
        elements = _parse_formula(formula)
        if not elements:
            raise ValueError("Не удалось разобрать формулу")

        # Определяем степени окисления
        oxidation = {}
        remaining = dict(elements)

        # Кислород обычно -2
        if "O" in remaining:
            oxidation["O"] = -2
            del remaining["O"]

        # Водород обычно +1 (кроме гидридов металлов)
        if "H" in remaining:
            oxidation["H"] = 1
            del remaining["H"]

        # Щелочные металлы +1
        for elem in list(remaining.keys()):
            if elem in ALKALI_METALS:
                oxidation[elem] = 1
                del remaining[elem]

        # Щёлочноземельные металлы +2
        for elem in list(remaining.keys()):
            if elem in ALKALINE_EARTH:
                oxidation[elem] = 2
                del remaining[elem]

        # Фтор всегда -1
        if "F" in remaining:
            oxidation["F"] = -1
            del remaining["F"]

        # Для оставшихся элементов вычисляем по балансу зарядов
        if remaining:
            # Сумма известных степеней окисления * количество
            known_sum = sum(oxidation.get(elem, 0) * elements.get(elem, 0) for elem in oxidation)
            # Оставшиеся элементы
            remaining_count = sum(remaining.values())
            if len(remaining) == 1:
                elem = list(remaining.keys())[0]
                count = remaining[elem]
                if count > 0:
                    oxidation[elem] = -known_sum // count if known_sum != 0 else 0
            else:
                # Для нескольких элементов используем типичные степени окисления
                for elem in remaining:
                    typical = OXIDATION_STATES.get(elem, [0])
                    oxidation[elem] = typical[0] if typical else 0

        return {
            "formula": formula,
            "oxidation_states": oxidation,
            "elements": elements,
        }
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


@router.get("/solubility", summary="Проверка растворимости")
def check_solubility(
    cation: str = Query(..., description="Катион (например, Na, Ca, Ag)"),
    anion: str = Query(..., description="Анион (например, Cl, SO4, CO3)"),
    lang: str = Depends(get_lang),
):
    """
    Проверяет растворимость соединения по таблице растворимости.
    Возвращает: soluble (растворим), insoluble (нерастворим), slightly (малорастворим).
    """
    try:
        # Нитраты, ацетаты, хлораты, перхлораты всегда растворимы
        if anion in ("NO3", "CH3COO", "ClO3", "ClO4"):
            return {
                "compound": f"{cation}{anion}",
                "solubility": "soluble",
                "solubility_ru": "растворим",
                "rule": "Все нитраты, ацетаты, хлораты и перхлораты растворимы",
            }

        # Хлориды, бромиды, иодиды
        if anion in ("Cl", "Br", "I"):
            exceptions = SOLUBILITY_TABLE.get(anion, {})
            if cation in exceptions:
                sol = exceptions[cation]
            else:
                sol = "soluble"
            return {
                "compound": f"{cation}{anion}",
                "solubility": sol,
                "solubility_ru": {"soluble": "растворим", "insoluble": "нерастворим", "slightly": "малорастворим"}.get(sol, sol),
                "rule": f"Исключение для {cation}" if cation in exceptions else "Все галогениды растворимы, кроме Ag, Pb, Hg2, Tl",
                "rule_en": f"Exception for {cation}" if cation in exceptions else "All halides soluble except Ag, Pb, Hg2, Tl",
                "rule_uk": f"Виняток для {cation}" if cation in exceptions else "Усі галогеніди розчинні, крім Ag, Pb, Hg2, Tl",
                "rule_es": f"Excepción para {cation}" if cation in exceptions else "Todos los haluros solubles excepto Ag, Pb, Hg2, Tl",
            }

        # Сульфаты
        if anion == "SO4":
            exceptions = SOLUBILITY_TABLE.get("SO4", {})
            if cation in exceptions:
                sol = exceptions[cation]
            else:
                sol = "soluble"
            return {
                "compound": f"{cation}SO4",
                "solubility": sol,
                "solubility_ru": {"soluble": "растворим", "insoluble": "нерастворим", "slightly": "малорастворим"}.get(sol, sol),
                "rule": f"Исключение для {cation}" if cation in exceptions else "Все сульфаты растворимы, кроме Ba, Pb, Sr, Ca(мал)",
                "rule_en": f"Exception for {cation}" if cation in exceptions else "All sulfates soluble except Ba, Pb, Sr, Ca(slightly)",
                "rule_uk": f"Виняток для {cation}" if cation in exceptions else "Усі сульфати розчинні, крім Ba, Pb, Sr, Ca(мала)",
                "rule_es": f"Excepción para {cation}" if cation in exceptions else "Todos los sulfatos solubles excepto Ba, Pb, Sr, Ca(poco)",
            }

        # Сульфиды
        if anion == "S":
            if cation in ALKALI_METALS or cation == "NH4":
                sol = "soluble"
            elif cation in ALKALINE_EARTH:
                sol = "slightly"
            else:
                sol = "insoluble"
            return {
                "compound": f"{cation}S",
                "solubility": sol,
                "solubility_ru": {"soluble": "растворим", "insoluble": "нерастворим", "slightly": "малорастворим"}.get(sol, sol),
                "rule": "Сульфиды щелочных и щёлочноземельных металлов растворимы",
                "rule_en": "Sulfides of alkali and alkaline earth metals are soluble",
                "rule_uk": "Сульфіди лужних та лужноземельних металів розчинні",
                "rule_es": "Sulfuros de metales alcalinos y alcalinotérreos son solubles",
            }

        # Карбонаты
        if anion == "CO3":
            if cation in ALKALI_METALS or cation == "NH4":
                sol = "soluble"
            elif cation == "Mg":
                sol = "slightly"
            else:
                sol = "insoluble"
            return {
                "compound": f"{cation}CO3",
                "solubility": sol,
                "solubility_ru": {"soluble": "растворим", "insoluble": "нерастворим", "slightly": "малорастворим"}.get(sol, sol),
                "rule": "Карбонаты щелочных металлов растворимы, остальные — нет",
                "rule_en": "Carbonates of alkali metals are soluble, others are not",
                "rule_uk": "Карбонати лужних металів розчинні, інші — ні",
                "rule_es": "Carbonatos de metales alcalinos son solubles, otros no",
            }

        # Гидроксиды
        if anion == "OH":
            if cation in ALKALI_METALS or cation == "NH4":
                sol = "soluble"
            elif cation in ("Ba", "Sr"):
                sol = "soluble"
            elif cation == "Ca":
                sol = "slightly"
            else:
                sol = "insoluble"
            return {
                "compound": f"{cation}OH",
                "solubility": sol,
                "solubility_ru": {"soluble": "растворим", "insoluble": "нерастворим", "slightly": "малорастворим"}.get(sol, sol),
                "rule": "Гидроксиды щелочных металлов, Ba, Sr растворимы",
                "rule_en": "Hydroxides of alkali metals, Ba, Sr are soluble",
                "rule_uk": "Гідроксиди лужних металів, Ba, Sr розчинні",
                "rule_es": "Hidróxidos de metales alcalinos, Ba, Sr son solubles",
            }

        # Фосфаты
        if anion == "PO4":
            if cation in ALKALI_METALS or cation == "NH4":
                sol = "soluble"
            elif cation == "Mg":
                sol = "slightly"
            else:
                sol = "insoluble"
            return {
                "compound": f"{cation}PO4",
                "solubility": sol,
                "solubility_ru": {"soluble": "растворим", "insoluble": "нерастворим", "slightly": "малорастворим"}.get(sol, sol),
                "rule": "Фосфаты щелочных металлов растворимы",
                "rule_en": "Phosphates of alkali metals are soluble",
                "rule_uk": "Фосфати лужних металів розчинні",
                "rule_es": "Fosfatos de metales alcalinos son solubles",
            }

        # По умолчанию
        return {
            "compound": f"{cation}{anion}",
            "solubility": "unknown",
            "solubility_ru": "неизвестно",
            "rule": "Нет данных в таблице растворимости",
        }
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


@router.get("/redox", summary="Окислительно-восстановительные реакции")
def analyze_redox(
    reaction: str = Query(..., description="Реакция в формате: Fe2O3 + CO -> Fe + CO2"),
    lang: str = Depends(get_lang),
):
    """
    Анализирует ОВР:
    - Определяет окислитель и восстановитель
    - Уравнивает реакцию
    - Показывает перенос электронов
    """
    try:
        # Парсим реакцию
        parts = reaction.split("->")
        if len(parts) != 2:
            raise ValueError("Реакция должна быть в формате: A + B -> C + D")

        left = [s.strip() for s in parts[0].split("+")]
        right = [s.strip() for s in parts[1].split("+")]

        # Парсим формулы
        left_formulas = [_parse_formula(s) for s in left]
        right_formulas = [_parse_formula(s) for s in right]

        # Собираем все элементы
        all_elements = set()
        for f in left_formulas + right_formulas:
            all_elements.update(f.keys())

        # Уравниваем реакцию
        balanced = _balance_reaction(left, right, left_formulas, right_formulas, all_elements)

        return {
            "reaction": reaction,
            "balanced": balanced["equation"],
            "coefficients": balanced["coefficients"],
            "oxidizing_agent": balanced.get("oxidizing_agent"),
            "reducing_agent": balanced.get("reducing_agent"),
            "electron_transfer": balanced.get("electron_transfer"),
        }
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


def _balance_reaction(left_formulas: List[str], right_formulas: List[str], left_parsed: List[Dict], right_parsed: List[Dict], elements: set) -> dict:
    """
    Уравнивает химическую реакцию перебором минимальных целых коэффициентов.
    """
    from math import gcd
    from functools import reduce
    from itertools import product

    n_left = len(left_formulas)
    n_right = len(right_formulas)
    n_compounds = n_left + n_right
    elem_list = sorted(elements)

    if n_compounds < 2:
        return {"equation": " + ".join(["?"] * n_compounds), "coefficients": []}

    # Перебор коэффициентов от 1 до 20
    max_coeff = 20
    for coeffs in product(range(1, max_coeff + 1), repeat=n_compounds):
        # Проверяем баланс по каждому элементу
        balanced = True
        for elem in elem_list:
            left_count = sum(coeffs[j] * left_parsed[j].get(elem, 0) for j in range(n_left))
            right_count = sum(coeffs[n_left + j] * right_parsed[j].get(elem, 0) for j in range(n_right))
            if left_count != right_count:
                balanced = False
                break

        if balanced:
            # Сокращаем на НОД
            common = reduce(gcd, coeffs)
            int_coeffs = [c // common for c in coeffs]

            # Формируем уравнение
            left_str = " + ".join(f"{int_coeffs[i]}{left_formulas[i]}" if int_coeffs[i] != 1 else left_formulas[i] for i in range(n_left))
            right_str = " + ".join(f"{int_coeffs[n_left + i]}{right_formulas[i]}" if int_coeffs[n_left + i] != 1 else right_formulas[i] for i in range(n_right))
            equation = f"{left_str} -> {right_str}"

            return {
                "equation": equation,
                "coefficients": int_coeffs,
            }

    return {"equation": "Не удалось уравнять", "coefficients": [], "equation_en": "Could not balance", "equation_uk": "Не вдалося врівняти", "equation_es": "No se pudo balancear"}


@router.get("/stoichiometry", summary="Стехиометрия реакции")
def calculate_stoichiometry(
    reaction: str = Query(..., description="Реакция: 2H2 + O2 -> 2H2O"),
    given_mass: Optional[float] = Query(None, description="Масса реагента в г"),
    given_moles: Optional[float] = Query(None, description="Количество молей реагента"),
    lang: str = Depends(get_lang),
):
    """
    Рассчитывает стехиометрию реакции:
    - Массы всех продуктов
    - Количество молей всех веществ
    - Теоретический выход
    """
    try:
        parts = reaction.split("->")
        if len(parts) != 2:
            raise ValueError("Реакция должна быть в формате: A + B -> C + D")

        left = [s.strip() for s in parts[0].split("+")]
        right = [s.strip() for s in parts[1].split("+")]

        # Парсим формулы с коэффициентами
        left_parsed = []
        for s in left:
            # Проверяем коэффициент
            coeff = 1
            formula = s
            if s[0].isdigit():
                i = 0
                while i < len(s) and s[i].isdigit():
                    i += 1
                coeff = int(s[:i])
                formula = s[i:]
            left_parsed.append({"formula": formula, "coeff": coeff, "elements": _parse_formula(formula)})

        right_parsed = []
        for s in right:
            coeff = 1
            formula = s
            if s[0].isdigit():
                i = 0
                while i < len(s) and s[i].isdigit():
                    i += 1
                coeff = int(s[:i])
                formula = s[i:]
            right_parsed.append({"formula": formula, "coeff": coeff, "elements": _parse_formula(formula)})

        # Вычисляем молекулярные массы
        for item in left_parsed + right_parsed:
            item["molar_mass"] = sum(ATOMIC_MASSES.get(e, 0) * c for e, c in item["elements"].items())

        # Если дана масса или моли реагента
        if given_mass is not None or given_moles is not None:
            # Берём первый реагент
            reactant = left_parsed[0]
            if given_moles is None:
                given_moles = given_mass / reactant["molar_mass"]

            # Мольное соотношение
            ratio = given_moles / reactant["coeff"]

            # Рассчитываем продукты
            products = []
            for p in right_parsed:
                moles = ratio * p["coeff"]
                mass = moles * p["molar_mass"]
                products.append({
                    "formula": p["formula"],
                    "moles": round(moles, 4),
                    "mass_g": round(mass, 4),
                })

            return {
                "reaction": reaction,
                "given": {
                    "reactant": reactant["formula"],
                    "mass_g": given_mass,
                    "moles": round(given_moles, 4),
                },
                "products": products,
            }

        # Просто возвращаем стехиометрию
        return {
            "reaction": reaction,
            "reactants": [{"formula": r["formula"], "coeff": r["coeff"], "molar_mass": round(r["molar_mass"], 4)} for r in left_parsed],
            "products": [{"formula": p["formula"], "coeff": p["coeff"], "molar_mass": round(p["molar_mass"], 4)} for p in right_parsed],
        }
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


@router.get("/orbitals", summary="Электронная конфигурация и орбитали")
def get_orbitals(
    element: str = Query(..., description="Символ элемента (например, Fe)"),
    lang: str = Depends(get_lang),
):
    """
    Возвращает электронную конфигурацию элемента:
    - Конфигурацию в нотации благородного газа
    - Распределение по орбиталям
    - Число неспаренных электронов
    - Валентные электроны
    """
    try:
        # Атомные номера
        atomic_numbers = {
            "H": 1, "He": 2, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9, "Ne": 10,
            "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15, "S": 16, "Cl": 17, "Ar": 18,
            "K": 19, "Ca": 20, "Sc": 21, "Ti": 22, "V": 23, "Cr": 24, "Mn": 25, "Fe": 26,
            "Co": 27, "Ni": 28, "Cu": 29, "Zn": 30, "Ga": 31, "Ge": 32, "As": 33, "Se": 34,
            "Br": 35, "Kr": 36, "Rb": 37, "Sr": 38, "Y": 39, "Zr": 40, "Nb": 41, "Mo": 42,
            "Tc": 43, "Ru": 44, "Rh": 45, "Pd": 46, "Ag": 47, "Cd": 48, "In": 49, "Sn": 50,
            "Sb": 51, "Te": 52, "I": 53, "Xe": 54, "Cs": 55, "Ba": 56, "La": 57, "Ce": 58,
            "Pr": 59, "Nd": 60, "Pm": 61, "Sm": 62, "Eu": 63, "Gd": 64, "Tb": 65, "Dy": 66,
            "Ho": 67, "Er": 68, "Tm": 69, "Yb": 70, "Lu": 71, "Hf": 72, "Ta": 73, "W": 74,
            "Re": 75, "Os": 76, "Ir": 77, "Pt": 78, "Au": 79, "Hg": 80, "Tl": 81, "Pb": 82,
            "Bi": 83, "Po": 84, "At": 85, "Rn": 86, "Fr": 87, "Ra": 88, "Ac": 89, "Th": 90,
            "Pa": 91, "U": 92, "Np": 93, "Pu": 94, "Am": 95, "Cm": 96, "Bk": 97, "Cf": 98,
            "Es": 99, "Fm": 100, "Md": 101, "No": 102, "Lr": 103, "Rf": 104, "Db": 105,
            "Sg": 106, "Bh": 107, "Hs": 108, "Mt": 109, "Ds": 110, "Rg": 111, "Cn": 112,
            "Nh": 113, "Fl": 114, "Mc": 115, "Lv": 116, "Ts": 117, "Og": 118,
        }

        Z = atomic_numbers.get(element)
        if Z is None:
            raise ValueError(f"Неизвестный элемент: {element}")

        # Порядок заполнения орбиталей
        orbital_order = [
            ("1s", 2), ("2s", 2), ("2p", 6), ("3s", 2), ("3p", 6), ("4s", 2), ("3d", 10),
            ("4p", 6), ("5s", 2), ("4d", 10), ("5p", 6), ("6s", 2), ("4f", 14), ("5d", 10),
            ("6p", 6), ("7s", 2), ("5f", 14), ("6d", 10), ("7p", 6),
        ]

        # Заполняем орбитали
        remaining = Z
        config = []
        for orbital, capacity in orbital_order:
            if remaining <= 0:
                break
            electrons = min(remaining, capacity)
            config.append(f"{orbital}{electrons}")
            remaining -= electrons

        # Нотация благородного газа
        noble_gases = {2: "He", 10: "Ne", 18: "Ar", 36: "Kr", 54: "Xe", 86: "Rn", 118: "Og"}
        core = ""
        core_Z = 0
        for ng_Z, ng in sorted(noble_gases.items()):
            if ng_Z < Z:
                core = ng
                core_Z = ng_Z

        if core:
            # Упрощённо: берём всё после конфигурации благородного газа
            ng_config = _get_noble_gas_config(core)
            ng_orbitals = len(ng_config.split())
            valence_config = " ".join(config[ng_orbitals:])
            config_str = f"[{core}] {valence_config}"
        else:
            config_str = " ".join(config)

        # Неспарённые электроны (по правилу Хунда)
        unpaired = _count_unpaired(config)

        # Валентные электроны
        valence = _count_valence(config)

        return {
            "element": element,
            "atomic_number": Z,
            "electron_configuration": config_str,
            "full_configuration": " ".join(config),
            "unpaired_electrons": unpaired,
            "valence_electrons": valence,
            "orbitals": [{"orbital": o, "electrons": int(c)} for o, c in [c.split(":") if ":" in c else (c[:-1], c[-1]) for c in config]],
        }
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


def _get_noble_gas_config(element: str) -> str:
    """Возвращает конфигурацию благородного газа."""
    configs = {
        "He": "1s2",
        "Ne": "1s2 2s2 2p6",
        "Ar": "1s2 2s2 2p6 3s2 3p6",
        "Kr": "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6",
        "Xe": "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6 5s2 4d10 5p6",
        "Rn": "1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6 5s2 4d10 5p6 6s2 4f14 5d10 6p6",
    }
    return configs.get(element, "")


def _count_unpaired(config: list) -> int:
    """Подсчитывает неспаренные электроны по правилу Хунда."""
    unpaired = 0
    for item in config:
        orbital = item[:-1]
        electrons = int(item[-1])
        if orbital.endswith("s"):
            unpaired += 0 if electrons == 2 else 1
        elif orbital.endswith("p"):
            # p-орбиталь: 3 подорбитали, каждая до 2 электронов
            if electrons <= 3:
                unpaired += electrons
            else:
                unpaired += 6 - electrons
        elif orbital.endswith("d"):
            # d-орбиталь: 5 подорбиталей
            if electrons <= 5:
                unpaired += electrons
            else:
                unpaired += 10 - electrons
        elif orbital.endswith("f"):
            # f-орбиталь: 7 подорбиталей
            if electrons <= 7:
                unpaired += electrons
            else:
                unpaired += 14 - electrons
    return unpaired


def _count_valence(config: list) -> int:
    """Подсчитывает валентные электроны (последняя оболочка)."""
    if not config:
        return 0
    # Последняя оболочка
    last_n = int(config[-1][0])
    valence = 0
    for item in config:
        n = int(item[0])
        if n == last_n:
            valence += int(item[-1])
        # Также учитываем d-электронны предыдущей оболочки для переходных металлов
        if n == last_n - 1 and item[1] == "d":
            valence += int(item[-1])
    return valence


@router.get("/element_info", summary="Информация об элементе")
def get_element_info(
    symbol: str = Query(..., description="Символ элемента (например, Fe)"),
    lang: str = Depends(get_lang),
):
    """
    Получает информацию об элементе:
    - Атомный номер и масса
    - Электроотрицательность
    - Степени окисления
    - Группа и период
    """
    try:
        # Получаем из PubChem
        info = _get_element_info(symbol)

        # Добавляем локальные данные
        info["electronegativity"] = ELECTRONEGATIVITY.get(symbol)
        info["oxidation_states"] = OXIDATION_STATES.get(symbol, [])
        info["atomic_mass"] = ATOMIC_MASSES.get(symbol)

        # Группа и период (упрощённо)
        atomic_numbers = {
            "H": 1, "He": 2, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9, "Ne": 10,
            "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15, "S": 16, "Cl": 17, "Ar": 18,
            "K": 19, "Ca": 20, "Sc": 21, "Ti": 22, "V": 23, "Cr": 24, "Mn": 25, "Fe": 26,
            "Co": 27, "Ni": 28, "Cu": 29, "Zn": 30, "Ga": 31, "Ge": 32, "As": 33, "Se": 34,
            "Br": 35, "Kr": 36, "Rb": 37, "Sr": 38, "Y": 39, "Zr": 40, "Nb": 41, "Mo": 42,
            "Tc": 43, "Ru": 44, "Rh": 45, "Pd": 46, "Ag": 47, "Cd": 48, "In": 49, "Sn": 50,
            "Sb": 51, "Te": 52, "I": 53, "Xe": 54, "Cs": 55, "Ba": 56, "La": 57, "Ce": 58,
            "Pr": 59, "Nd": 60, "Pm": 61, "Sm": 62, "Eu": 63, "Gd": 64, "Tb": 65, "Dy": 66,
            "Ho": 67, "Er": 68, "Tm": 69, "Yb": 70, "Lu": 71, "Hf": 72, "Ta": 73, "W": 74,
            "Re": 75, "Os": 76, "Ir": 77, "Pt": 78, "Au": 79, "Hg": 80, "Tl": 81, "Pb": 82,
            "Bi": 83, "Po": 84, "At": 85, "Rn": 86, "Fr": 87, "Ra": 88, "Ac": 89, "Th": 90,
            "Pa": 91, "U": 92, "Np": 93, "Pu": 94, "Am": 95, "Cm": 96, "Bk": 97, "Cf": 98,
            "Es": 99, "Fm": 100, "Md": 101, "No": 102, "Lr": 103, "Rf": 104, "Db": 105,
            "Sg": 106, "Bh": 107, "Hs": 108, "Mt": 109, "Ds": 110, "Rg": 111, "Cn": 112,
            "Nh": 113, "Fl": 114, "Mc": 115, "Lv": 116, "Ts": 117, "Og": 118,
        }
        Z = atomic_numbers.get(symbol)
        if Z:
            info["atomic_number"] = Z
            # Период
            if Z <= 2: info["period"] = 1
            elif Z <= 10: info["period"] = 2
            elif Z <= 18: info["period"] = 3
            elif Z <= 36: info["period"] = 4
            elif Z <= 54: info["period"] = 5
            elif Z <= 86: info["period"] = 6
            else: info["period"] = 7

        return info
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))
