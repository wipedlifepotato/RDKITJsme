"""
Интернационализация (i18n) для SMILES API Rdkit.
Поддерживаемые языки: ru (по умолчанию), uk, en, es.
"""

from typing import Dict, Optional
from fastapi import Query

# Поддерживаемые языки
SUPPORTED_LANGUAGES = {"ru", "uk", "en", "es"}
DEFAULT_LANGUAGE = "ru"

# Переводы
TRANSLATIONS: Dict[str, Dict[str, str]] = {
    # ═══════════════════════════════════════════════════════════════════════════
    # Общие сообщения об ошибках
    # ═══════════════════════════════════════════════════════════════════════════
    "invalid_smiles": {
        "ru": "Некорректная строка SMILES",
        "uk": "Некоректний рядок SMILES",
        "en": "Invalid SMILES string",
        "es": "Cadena SMILES inválida",
    },
    "invalid_smiles_structure": {
        "ru": "Некорректная структура SMILES",
        "uk": "Некоректна структура SMILES",
        "en": "Invalid SMILES structure",
        "es": "Estructura SMILES inválida",
    },
    "empty_smiles": {
        "ru": "Пустая строка SMILES",
        "uk": "Порожній рядок SMILES",
        "en": "Empty SMILES string",
        "es": "Cadena SMILES vacía",
    },
    "name_not_found": {
        "ru": "Название не найдено",
        "uk": "Назву не знайдено",
        "en": "Name not found",
        "es": "Nombre no encontrado",
    },
    "proxy_request_error": {
        "ru": "Ошибка прокси/запроса",
        "uk": "Помилка проксі/запиту",
        "en": "Proxy/Request error",
        "es": "Error de proxy/solicitud",
    },
    "molecule_not_found_in_library": {
        "ru": "Молекула не найдена в библиотеке",
        "uk": "Молекулу не знайдено в бібліотеці",
        "en": "Molecule not found in library",
        "es": "Molécula no encontrada en la biblioteca",
    },
    "unsupported_format": {
        "ru": "Неподдерживаемый формат",
        "uk": "Непідтримуваний формат",
        "en": "Unsupported format",
        "es": "Formato no soportado",
    },
    "failed_to_embed_molecule": {
        "ru": "Не удалось встроить молекулу в 3D пространство",
        "uk": "Не вдалося вбудувати молекулу в 3D простір",
        "en": "Failed to embed molecule in 3D space",
        "es": "No se pudo incrustar la molécula en el espacio 3D",
    },
    "failed_to_generate_conformers": {
        "ru": "Не удалось сгенерировать конформеры",
        "uk": "Не вдалося згенерувати конформери",
        "en": "Failed to generate conformers",
        "es": "No se pudieron generar confórmeros",
    },
    "need_at_least_two_smiles": {
        "ru": "Для сравнения нужно минимум 2 корректных SMILES",
        "uk": "Для порівняння потрібно мінімум 2 коректних SMILES",
        "en": "Need at least 2 valid SMILES for comparison",
        "es": "Se necesitan al menos 2 SMILES válidos para comparar",
    },
    "invalid_inchikey_format": {
        "ru": "Некорректный формат InChIKey. Ожидаемый формат: XXXXXXXXXXXXXX-XXXXXXXXXX-X",
        "uk": "Некоректний формат InChIKey. Очікуваний формат: XXXXXXXXXXXXXX-XXXXXXXXXX-X",
        "en": "Invalid InChIKey format. Expected format: XXXXXXXXXXXXXX-XXXXXXXXXX-X",
        "es": "Formato InChIKey inválido. Formato esperado: XXXXXXXXXXXXXX-XXXXXXXXXX-X",
    },
    "could_not_resolve_inchikey": {
        "ru": "Не удалось разрешить InChIKey",
        "uk": "Не вдалося розв'язати InChIKey",
        "en": "Could not resolve InChIKey",
        "es": "No se pudo resolver InChIKey",
    },
    "resolved_smiles_invalid": {
        "ru": "Разрешённый SMILES некорректен",
        "uk": "Розв'язаний SMILES некоректний",
        "en": "Resolved SMILES is invalid",
        "es": "El SMILES resuelto es inválido",
    },
    "invalid_smarts_pattern": {
        "ru": "Некорректный SMARTS-паттерн",
        "uk": "Некоректний SMARTS-шаблон",
        "en": "Invalid SMARTS pattern",
        "es": "Patrón SMARTS inválido",
    },
    "invalid_smirks_reaction": {
        "ru": "Некорректный SMIRKS-паттерн реакции",
        "uk": "Некоректний SMIRKS-шаблон реакції",
        "en": "Invalid SMIRKS reaction pattern",
        "es": "Patrón de reacción SMIRKS inválido",
    },
    "no_reaction_products": {
        "ru": "Продукты реакции не образованы",
        "uk": "Продукти реакції не утворилися",
        "en": "No reaction products formed",
        "es": "No se formaron productos de reacción",
    },
    "invalid_smiles_or_smarts": {
        "ru": "Некорректные SMILES или SMARTS-паттерн",
        "uk": "Некоректні SMILES або SMARTS-шаблон",
        "en": "Invalid SMILES or SMARTS pattern",
        "es": "SMILES o patrón SMARTS inválido",
    },

    # ═══════════════════════════════════════════════════════════════════════════
    # Описания эндпоинтов
    # ═══════════════════════════════════════════════════════════════════════════
    "ep_similarity_desc": {
        "ru": "Вычисляет сходство Танимото между двумя молекулами",
        "uk": "Обчислює схожість Танімото між двома молекулами",
        "en": "Calculate Tanimoto similarity between two molecules",
        "es": "Calcular la similitud de Tanimoto entre dos moléculas",
    },
    "ep_3d_sdf_desc": {
        "ru": "Генерирует 3D-структуру молекулы в формате SDF",
        "uk": "Генерує 3D-структуру молекули у форматі SDF",
        "en": "Generate 3D structure of molecule in SDF format",
        "es": "Generar estructura 3D de la molécula en formato SDF",
    },
    "ep_get_name_desc": {
        "ru": "Получает IUPAC-название молекулы по SMILES",
        "uk": "Отримує IUPAC-назву молекули за SMILES",
        "en": "Get IUPAC name of molecule by SMILES",
        "es": "Obtener el nombre IUPAC de la molécula por SMILES",
    },
    "ep_get_chiral_desc": {
        "ru": "Находит хиральные центры в молекуле",
        "uk": "Знаходить хіральні центри в молекулі",
        "en": "Find chiral centers in molecule",
        "es": "Encontrar centros quirales en la molécula",
    },
    "ep_get_properties_desc": {
        "ru": "Вычисляет молекулярные свойства (формула, масса, LogP и др.)",
        "uk": "Обчислює молекулярні властивості (формула, маса, LogP тощо)",
        "en": "Calculate molecular properties (formula, mass, LogP, etc.)",
        "es": "Calcular propiedades moleculares (fórmula, masa, LogP, etc.)",
    },
    "ep_convert_desc": {
        "ru": "Конвертирует SMILES в канонический, изомерный, InChI и InChIKey",
        "uk": "Конвертує SMILES у канонічний, ізомерний, InChI та InChIKey",
        "en": "Convert SMILES to canonical, isomeric, InChI and InChIKey",
        "es": "Convertir SMILES a canónico, isomérico, InChI e InChIKey",
    },
    "ep_substructure_desc": {
        "ru": "Поиск подструктуры по SMARTS-паттерну в молекуле",
        "uk": "Пошук підструктури за SMARTS-шаблоном у молекулі",
        "en": "Search substructure by SMARTS pattern in molecule",
        "es": "Buscar subestructura por patrón SMARTS en la molécula",
    },
    "ep_render_desc": {
        "ru": "Рендерит 2D-изображение молекулы",
        "uk": "Рендерить 2D-зображення молекули",
        "en": "Render 2D image of molecule",
        "es": "Renderizar imagen 2D de la molécula",
    },
    "ep_adme_desc": {
        "ru": "ADME-свойства, Липинский, QED, SA Score",
        "uk": "ADME-властивості, Ліпінський, QED, SA Score",
        "en": "ADME properties, Lipinski, QED, SA Score",
        "es": "Propiedades ADME, Lipinski, QED, SA Score",
    },
    "ep_functional_groups_desc": {
        "ru": "Поиск функциональных групп (25+ типов)",
        "uk": "Пошук функціональних груп (25+ типів)",
        "en": "Search functional groups (25+ types)",
        "es": "Buscar grupos funcionales (25+ tipos)",
    },
    "ep_validate_desc": {
        "ru": "Валидация SMILES с диагностикой",
        "uk": "Валідація SMILES з діагностикою",
        "en": "SMILES validation with diagnostics",
        "es": "Validación SMILES con diagnóstico",
    },
    "ep_pka_desc": {
        "ru": "Оценка pKa (кислотные/основные центры)",
        "uk": "Оцінка pKa (кислотні/основні центри)",
        "en": "Estimate pKa (acidic/basic centers)",
        "es": "Estimar pKa (centros ácidos/básicos)",
    },
    "ep_descriptors_desc": {
        "ru": "Все дескрипторы RDKit (200+)",
        "uk": "Всі дескриптори RDKit (200+)",
        "en": "All RDKit descriptors (200+)",
        "es": "Todos los descriptores RDKit (200+)",
    },
    "ep_reaction_desc": {
        "ru": "Применение реакций (SMIRKS)",
        "uk": "Застосування реакцій (SMIRKS)",
        "en": "Apply reactions (SMIRKS)",
        "es": "Aplicar reacciones (SMIRKS)",
    },
    "ep_by_inchikey_desc": {
        "ru": "Поиск по InChIKey",
        "uk": "Пошук за InChIKey",
        "en": "Search by InChIKey",
        "es": "Buscar por InChIKey",
    },
    "ep_3d_desc": {
        "ru": "3D-структура (SDF/MOL/PDB)",
        "uk": "3D-структура (SDF/MOL/PDB)",
        "en": "3D structure (SDF/MOL/PDB)",
        "es": "Estructura 3D (SDF/MOL/PDB)",
    },
    "ep_render_highlighted_desc": {
        "ru": "Рендеринг с подсветкой групп",
        "uk": "Рендеринг з підсвіткою груп",
        "en": "Render with group highlighting",
        "es": "Renderizar con resaltado de grupos",
    },
    "ep_batch_desc": {
        "ru": "Пакетный анализ молекул",
        "uk": "Пакетний аналіз молекул",
        "en": "Batch molecule analysis",
        "es": "Análisis de moléculas por lotes",
    },
    "ep_library_desc": {
        "ru": "Список молекул в библиотеке",
        "uk": "Список молекул у бібліотеці",
        "en": "List of molecules in library",
        "es": "Lista de moléculas en la biblioteca",
    },
    "ep_library_save_desc": {
        "ru": "Сохранение в библиотеку",
        "uk": "Збереження в бібліотеку",
        "en": "Save to library",
        "es": "Guardar en biblioteca",
    },
    "ep_library_delete_desc": {
        "ru": "Удаление из библиотеки",
        "uk": "Видалення з бібліотеки",
        "en": "Delete from library",
        "es": "Eliminar de la biblioteca",
    },
    "ep_compare_desc": {
        "ru": "Сравнение нескольких молекул (матрица Tanimoto)",
        "uk": "Порівняння кількох молекул (матриця Tanimoto)",
        "en": "Compare multiple molecules (Tanimoto matrix)",
        "es": "Comparar múltiples moléculas (matriz Tanimoto)",
    },
    "ep_nmr_desc": {
        "ru": "Оценка химических сдвигов ЯМР (1H и 13C)",
        "uk": "Оцінка хімічних зсувів ЯМР (1H та 13C)",
        "en": "Estimate NMR chemical shifts (1H and 13C)",
        "es": "Estimar desplazamientos químicos RMN (1H y 13C)",
    },
    "ep_toxicity_desc": {
        "ru": "Оценка токсичности через каталоги RDKit (Brenk / PAINS)",
        "uk": "Оцінка токсичності через каталоги RDKit (Brenk / PAINS)",
        "en": "Assess toxicity via RDKit catalogs (Brenk / PAINS)",
        "es": "Evaluar toxicidad mediante catálogos RDKit (Brenk / PAINS)",
    },
    "ep_retrosynthesis_desc": {
        "ru": "Ретросинтетический анализ (RECAP)",
        "uk": "Ретросинтетичний аналіз (RECAP)",
        "en": "Retrosynthetic analysis (RECAP)",
        "es": "Análisis retrosintético (RECAP)",
    },
    "ep_conformers_desc": {
        "ru": "Генерация конформеров молекулы",
        "uk": "Генерація конформерів молекули",
        "en": "Generate molecule conformers",
        "es": "Generar confórmeros de la molécula",
    },
    "ep_export_svg_desc": {
        "ru": "Экспорт 2D-структуры в SVG",
        "uk": "Експорт 2D-структури в SVG",
        "en": "Export 2D structure to SVG",
        "es": "Exportar estructura 2D a SVG",
    },
    "ep_ms_desc": {
        "ru": "Масс-спектрометрия (изотопное распределение)",
        "uk": "Мас-спектрометрія (ізотопний розподіл)",
        "en": "Mass spectrometry (isotope distribution)",
        "es": "Espectrometría de masas (distribución isotópica)",
    },

    # ═══════════════════════════════════════════════════════════════════════════
    # Описания параметров
    # ═══════════════════════════════════════════════════════════════════════════
    "param_smiles": {
        "ru": "SMILES молекулы",
        "uk": "SMILES молекули",
        "en": "SMILES of molecule",
        "es": "SMILES de la molécula",
    },
    "param_smiles1": {
        "ru": "Первый SMILES",
        "uk": "Перший SMILES",
        "en": "First SMILES",
        "es": "Primer SMILES",
    },
    "param_smiles2": {
        "ru": "Второй SMILES",
        "uk": "Другий SMILES",
        "en": "Second SMILES",
        "es": "Segundo SMILES",
    },
    "param_smiles_to_name": {
        "ru": "SMILES для наименования",
        "uk": "SMILES для найменування",
        "en": "SMILES to name",
        "es": "SMILES para nombrar",
    },
    "param_smiles_to_convert": {
        "ru": "SMILES для конвертации",
        "uk": "SMILES для конвертації",
        "en": "SMILES to convert",
        "es": "SMILES para convertir",
    },
    "param_target_smiles": {
        "ru": "SMILES целевой молекулы",
        "uk": "SMILES цільової молекули",
        "en": "Target molecule SMILES",
        "es": "SMILES de la molécula objetivo",
    },
    "param_pattern_smarts": {
        "ru": "SMARTS-паттерн для поиска",
        "uk": "SMARTS-шаблон для пошуку",
        "en": "SMARTS pattern to search for",
        "es": "Patrón SMARTS para buscar",
    },
    "param_draw_type": {
        "ru": "Тип: just, indexes, charges, chiral",
        "uk": "Тип: just, indexes, charges, chiral",
        "en": "Type: just, indexes, charges, chiral",
        "es": "Tipo: just, indexes, charges, chiral",
    },
    "param_highlight_chiral": {
        "ru": "Подсветить хиральные центры",
        "uk": "Підсвітити хіральні центри",
        "en": "Highlight chiral centers",
        "es": "Resaltar centros quirales",
    },
    "param_reaction_smarts": {
        "ru": "SMARTS реакции (SMIRKS)",
        "uk": "SMARTS реакції (SMIRKS)",
        "en": "Reaction SMARTS (SMIRKS)",
        "es": "SMARTS de reacción (SMIRKS)",
    },
    "param_inchikey": {
        "ru": "InChIKey молекулы",
        "uk": "InChIKey молекули",
        "en": "InChIKey of molecule",
        "es": "InChIKey de la molécula",
    },
    "param_format": {
        "ru": "Формат: sdf, mol, pdb",
        "uk": "Формат: sdf, mol, pdb",
        "en": "Format: sdf, mol, pdb",
        "es": "Formato: sdf, mol, pdb",
    },
    "param_smiles_list": {
        "ru": "Список SMILES для сравнения",
        "uk": "Список SMILES для порівняння",
        "en": "List of SMILES for comparison",
        "es": "Lista de SMILES para comparación",
    },
    "param_highlight_groups": {
        "ru": "Ключи групп для подсветки",
        "uk": "Ключі груп для підсвітки",
        "en": "Group keys to highlight",
        "es": "Claves de grupos para resaltar",
    },
    "param_num_conformers": {
        "ru": "Количество конформеров",
        "uk": "Кількість конформерів",
        "en": "Number of conformers",
        "es": "Número de confórmeros",
    },
    "param_width": {
        "ru": "Ширина в пикселях",
        "uk": "Ширина у пікселях",
        "en": "Width in pixels",
        "es": "Ancho en píxeles",
    },
    "param_height": {
        "ru": "Высота в пикселях",
        "uk": "Висота у пікселях",
        "en": "Height in pixels",
        "es": "Alto en píxeles",
    },
    "param_smiles_to_validate": {
        "ru": "SMILES для проверки",
        "uk": "SMILES для перевірки",
        "en": "SMILES to validate",
        "es": "SMILES para validar",
    },
    "param_smiles_reagent": {
        "ru": "SMILES реагента",
        "uk": "SMILES реагента",
        "en": "SMILES of reagent",
        "es": "SMILES del reactivo",
    },
    "param_smiles_target": {
        "ru": "SMILES целевой молекулы",
        "uk": "SMILES цільової молекули",
        "en": "SMILES of target molecule",
        "es": "SMILES de la molécula objetivo",
    },

    # ═══════════════════════════════════════════════════════════════════════════
    # Названия функциональных групп
    # ═══════════════════════════════════════════════════════════════════════════
    "fg_hydroxyl": {
        "ru": "Гидроксильная (-OH)",
        "uk": "Гідроксильна (-OH)",
        "en": "Hydroxyl (-OH)",
        "es": "Hidroxilo (-OH)",
    },
    "fg_carboxylic_acid": {
        "ru": "Карбоновая кислота (-COOH)",
        "uk": "Карбонова кислота (-COOH)",
        "en": "Carboxylic acid (-COOH)",
        "es": "Ácido carboxílico (-COOH)",
    },
    "fg_amine_primary": {
        "ru": "Первичный амин (-NH2)",
        "uk": "Первинний амін (-NH2)",
        "en": "Primary amine (-NH2)",
        "es": "Amina primaria (-NH2)",
    },
    "fg_amine_secondary": {
        "ru": "Вторичный амин (-NH-)",
        "uk": "Вторинний амін (-NH-)",
        "en": "Secondary amine (-NH-)",
        "es": "Amina secundaria (-NH-)",
    },
    "fg_amine_tertiary": {
        "ru": "Третичный амин (-N<)",
        "uk": "Третинний амін (-N<)",
        "en": "Tertiary amine (-N<)",
        "es": "Amina terciaria (-N<)",
    },
    "fg_amide": {
        "ru": "Амид (-C(=O)N-)",
        "uk": "Амід (-C(=O)N-)",
        "en": "Amide (-C(=O)N-)",
        "es": "Amida (-C(=O)N-)",
    },
    "fg_nitro": {
        "ru": "Нитрогруппа (-NO2)",
        "uk": "Нітрупа (-NO2)",
        "en": "Nitro (-NO2)",
        "es": "Nitro (-NO2)",
    },
    "fg_sulfonyl": {
        "ru": "Сульфонильная (-SO2-)",
        "uk": "Сульфонільна (-SO2-)",
        "en": "Sulfonyl (-SO2-)",
        "es": "Sulfonilo (-SO2-)",
    },
    "fg_phosphate": {
        "ru": "Фосфат",
        "uk": "Фосфат",
        "en": "Phosphate",
        "es": "Fosfato",
    },
    "fg_aldehyde": {
        "ru": "Альдегидная (-CHO)",
        "uk": "Альдегідна (-CHO)",
        "en": "Aldehyde (-CHO)",
        "es": "Aldehído (-CHO)",
    },
    "fg_ketone": {
        "ru": "Кетоновая (>C=O)",
        "uk": "Кетонова (>C=O)",
        "en": "Ketone (>C=O)",
        "es": "Cetona (>C=O)",
    },
    "fg_ester": {
        "ru": "Сложноэфирная (-COOR)",
        "uk": "Сложноестерна (-COOR)",
        "en": "Ester (-COOR)",
        "es": "Éster (-COOR)",
    },
    "fg_ether": {
        "ru": "Простая эфирная (-O-)",
        "uk": "Проста ефірна (-O-)",
        "en": "Ether (-O-)",
        "es": "Éter (-O-)",
    },
    "fg_alkene": {
        "ru": "Алкеновая (C=C)",
        "uk": "Алкенова (C=C)",
        "en": "Alkene (C=C)",
        "es": "Alqueno (C=C)",
    },
    "fg_alkyne": {
        "ru": "Алкиновая (C≡C)",
        "uk": "Алкінова (C≡C)",
        "en": "Alkyne (C≡C)",
        "es": "Alquino (C≡C)",
    },
    "fg_aromatic_ring": {
        "ru": "Ароматическое кольцо (бензол)",
        "uk": "Ароматичне кільце (бензол)",
        "en": "Aromatic ring (benzene)",
        "es": "Anillo aromático (benceno)",
    },
    "fg_heterocycle": {
        "ru": "Гетероцикл",
        "uk": "Гетероцикл",
        "en": "Heterocycle",
        "es": "Heterociclo",
    },
    "fg_halogen_f": {
        "ru": "Фтор",
        "uk": "Фтор",
        "en": "Fluorine",
        "es": "Flúor",
    },
    "fg_halogen_cl": {
        "ru": "Хлор",
        "uk": "Хлор",
        "en": "Chlorine",
        "es": "Cloro",
    },
    "fg_halogen_br": {
        "ru": "Бром",
        "uk": "Бром",
        "en": "Bromine",
        "es": "Bromo",
    },
    "fg_halogen_i": {
        "ru": "Иод",
        "uk": "Йод",
        "en": "Iodine",
        "es": "Yodo",
    },
    "fg_thiol": {
        "ru": "Тиольная (-SH)",
        "uk": "Тіольна (-SH)",
        "en": "Thiol (-SH)",
        "es": "Tiol (-SH)",
    },
    "fg_nitrile": {
        "ru": "Нитрильная (-C≡N)",
        "uk": "Нітрильна (-C≡N)",
        "en": "Nitrile (-C≡N)",
        "es": "Nitrilo (-C≡N)",
    },
    "fg_azide": {
        "ru": "Азидная (-N3)",
        "uk": "Азидна (-N3)",
        "en": "Azide (-N3)",
        "es": "Azida (-N3)",
    },
    "fg_imine": {
        "ru": "Иминовая (C=N)",
        "uk": "Імінова (C=N)",
        "en": "Imine (C=N)",
        "es": "Imina (C=N)",
    },
    "fg_hydrazine": {
        "ru": "Гидразиновая (-N-N-)",
        "uk": "Гідразинова (-N-N-)",
        "en": "Hydrazine (-N-N-)",
        "es": "Hidracina (-N-N-)",
    },
    "fg_urea": {
        "ru": "Мочевина",
        "uk": "Сечовина",
        "en": "Urea",
        "es": "Urea",
    },
    "fg_carbamate": {
        "ru": "Карбамат",
        "uk": "Карбамат",
        "en": "Carbamate",
        "es": "Carbamato",
    },
    "fg_sulfonamide": {
        "ru": "Сульфонамид",
        "uk": "Сульфонамід",
        "en": "Sulfonamide",
        "es": "Sulfonamida",
    },

    # ═══════════════════════════════════════════════════════════════════════════
    # Валидация SMILES
    # ═══════════════════════════════════════════════════════════════════════════
    "val_empty": {
        "ru": "Строка SMILES пуста",
        "uk": "Рядок SMILES порожній",
        "en": "SMILES string is empty",
        "es": "La cadena SMILES está vacía",
    },
    "val_invalid_characters": {
        "ru": "Найдены недопустимые символы",
        "uk": "Знайдено неприпустимі символи",
        "en": "Invalid characters found",
        "es": "Caracteres inválidos encontrados",
    },
    "val_unbalanced_parentheses": {
        "ru": "Несбалансированные скобки",
        "uk": "Незбалансовані дужки",
        "en": "Unbalanced parentheses",
        "es": "Paréntesis desbalanceados",
    },
    "val_unbalanced_brackets": {
        "ru": "Несбалансированные квадратные скобки",
        "uk": "Незбалансовані квадратні дужки",
        "en": "Unbalanced brackets",
        "es": "Corchetes desbalanceados",
    },
    "val_valence_error": {
        "ru": "Ошибка валентности",
        "uk": "Помилка валентності",
        "en": "Valence error",
        "es": "Error de valencia",
    },
    "val_aromaticity_error": {
        "ru": "Ошибка ароматичности",
        "uk": "Помилка ароматичності",
        "en": "Aromaticity error",
        "es": "Error de aromaticidad",
    },
    "val_kekulization_error": {
        "ru": "Ошибка кекулизации",
        "uk": "Помилка кекулізації",
        "en": "Kekulization error",
        "es": "Error de kekulización",
    },
    "val_sanitization_error": {
        "ru": "Ошибка санитизации",
        "uk": "Помилка санітізації",
        "en": "Sanitization error",
        "es": "Error de sanitización",
    },
    "val_syntax_error": {
        "ru": "RDKit не смог разобрать эту строку SMILES",
        "uk": "RDKit не зміг розібрати цей рядок SMILES",
        "en": "RDKit could not parse this SMILES string",
        "es": "RDKit no pudo analizar esta cadena SMILES",
    },
    "val_large_molecule": {
        "ru": "Большая молекула (>100 атомов)",
        "uk": "Велика молекула (>100 атомів)",
        "en": "Large molecule (>100 atoms)",
        "es": "Molécula grande (>100 átomos)",
    },
    "val_high_mw": {
        "ru": "Высокая молекулярная масса (>1000 Да)",
        "uk": "Висока молекулярна маса (>1000 Да)",
        "en": "High molecular weight (>1000 Da)",
        "es": "Alto peso molecular (>1000 Da)",
    },
    "val_high_logp": {
        "ru": "Очень высокая липофильность (LogP > 7)",
        "uk": "Дуже висока пофільність (LogP > 7)",
        "en": "Very high lipophilicity (LogP > 7)",
        "es": "Muy alta lipofilicidad (LogP > 7)",
    },

    # ═══════════════════════════════════════════════════════════════════════════
    # QED интерпретация
    # ═══════════════════════════════════════════════════════════════════════════
    "qed_excellent": {
        "ru": "Отличный (drug-like)",
        "uk": "Відмінний (drug-like)",
        "en": "Excellent (drug-like)",
        "es": "Excelente (similar a fármaco)",
    },
    "qed_good": {
        "ru": "Хороший (умеренный drug-likeness)",
        "uk": "Гарний (помірний drug-likeness)",
        "en": "Good (moderate drug-likeness)",
        "es": "Bueno (semejanza moderada a fármaco)",
    },
    "qed_fair": {
        "ru": "Удовлетворительный (слабый drug-likeness)",
        "uk": "Задовільний (слабкий drug-likeness)",
        "en": "Fair (poor drug-likeness)",
        "es": "Regular (baja semejanza a fármaco)",
    },
    "qed_poor": {
        "ru": "Плохой (маловероятно drug-like)",
        "uk": "Поганий (малоймовірно drug-like)",
        "en": "Poor (unlikely drug-like)",
        "es": "Deficiente (poco probable de ser fármaco)",
    },

    # ═══════════════════════════════════════════════════════════════════════════
    # pKa типы центров
    # ═══════════════════════════════════════════════════════════════════════════
    "pka_carboxylic_acid": {
        "ru": "Карбоновая кислота",
        "uk": "Карбонова кислота",
        "en": "Carboxylic acid",
        "es": "Ácido carboxílico",
    },
    "pka_phenol": {
        "ru": "Фенол",
        "uk": "Фенол",
        "en": "Phenol",
        "es": "Fenol",
    },
    "pka_aliphatic_amine": {
        "ru": "Алифатический амин",
        "uk": "Аліфатичний амін",
        "en": "Aliphatic amine",
        "es": "Amina alifática",
    },
    "pka_aniline": {
        "ru": "Анилин",
        "uk": "Анілін",
        "en": "Aniline",
        "es": "Anilina",
    },
    "pka_thiol": {
        "ru": "Тиол",
        "uk": "Тіол",
        "en": "Thiol",
        "es": "Tiol",
    },
    "pka_sulfonic_acid": {
        "ru": "Сульфоновая кислота",
        "uk": "Сульфонова кислота",
        "en": "Sulfonic acid",
        "es": "Ácido sulfónico",
    },
    "pka_acidic": {
        "ru": "Кислотный",
        "uk": "Кислотний",
        "en": "Acidic",
        "es": "Ácido",
    },
    "pka_basic": {
        "ru": "Основной",
        "uk": "Основний",
        "en": "Basic",
        "es": "Básico",
    },
    "pka_strong_acid": {
        "ru": "Сильная кислота",
        "uk": "Сильна кислота",
        "en": "Strong acid",
        "es": "Ácido fuerte",
    },

    # ═══════════════════════════════════════════════════════════════════════════
    # Токсичность
    # ═══════════════════════════════════════════════════════════════════════════
    "tox_herg_name": {
        "ru": "Блокада hERG-канала",
        "uk": "Блокада hERG-каналу",
        "en": "hERG Channel Blockade",
        "es": "Bloqueo del canal hERG",
    },
    "tox_herg_desc": {
        "ru": "Блокирование калиевых каналов hERG — кардиотоксичность",
        "uk": "Блокування калієвих каналів hERG — кардіотоксичність",
        "en": "Blockade of hERG potassium channels — cardiotoxicity",
        "es": "Bloqueo de canales de potasio hERG — cardiotoxicidad",
    },
    "tox_herg_risk": {
        "ru": "Высокий риск удлинения интервала QT",
        "uk": "Високий ризик подовження інтервалу QT",
        "en": "High risk of QT interval prolongation",
        "es": "Alto riesgo de prolongación del intervalo QT",
    },
    "tox_ames_name": {
        "ru": "Тест Эймса (мутагенность)",
        "uk": "Тест Еймса (мутагенність)",
        "en": "Ames Test (Mutagenicity)",
        "es": "Test de Ames (Mutagenicidad)",
    },
    "tox_ames_desc": {
        "ru": "Мутагенность — повреждение ДНК",
        "uk": "Мутагенність — пошкодження ДНК",
        "en": "Mutagenicity — DNA damage",
        "es": "Mutagenicidad — daño al ADN",
    },
    "tox_ames_risk": {
        "ru": "Потенциальные мутагены",
        "uk": "Потенційні мутагени",
        "en": "Potential mutagens",
        "es": "Mutágenos potenciales",
    },
    "tox_hepato_name": {
        "ru": "Гепатотоксичность",
        "uk": "Гепатотоксичність",
        "en": "Hepatotoxicity",
        "es": "Hepatotoxicidad",
    },
    "tox_hepato_desc": {
        "ru": "Токсичность для печени",
        "uk": "Токсичність для печінки",
        "en": "Liver toxicity",
        "es": "Toxicidad hepática",
    },
    "tox_hepato_risk": {
        "ru": "Потенциальная гепатотоксичность",
        "uk": "Потенційна гепатотоксичність",
        "en": "Potential hepatotoxicity",
        "es": "Hepatotoxicidad potencial",
    },
    "tox_basic_amines": {
        "ru": "Основные амины",
        "uk": "Основні аміни",
        "en": "Basic amines",
        "es": "Aminas básicas",
    },
    "tox_anilines": {
        "ru": "Анилины",
        "uk": "Аніліни",
        "en": "Anilines",
        "es": "Anilinas",
    },
    "tox_pah": {
        "ru": "Полициклические ароматические углеводороды",
        "uk": "Поліциклічні ароматичні вуглеводні",
        "en": "Polycyclic aromatic hydrocarbons",
        "es": "Hidrocarburos aromáticos policíclicos",
    },
    "tox_azides": {
        "ru": "Азиды",
        "uk": "Азиди",
        "en": "Azides",
        "es": "Azidas",
    },
    "tox_nitro": {
        "ru": "Нитросоединения",
        "uk": "Нітросполуки",
        "en": "Nitro compounds",
        "es": "Compuestos nitro",
    },
    "tox_diamides": {
        "ru": "Диамиды",
        "uk": "Діаміди",
        "en": "Diamides",
        "es": "Diamidas",
    },
    "tox_sulfonic_acids": {
        "ru": "Сульфоновые кислоты",
        "uk": "Сульфонові кислоти",
        "en": "Sulfonic acids",
        "es": "Ácidos sulfónicos",
    },
    "tox_haloalkanes": {
        "ru": "Галогеналканы",
        "uk": "Галогеналкани",
        "en": "Haloalkanes",
        "es": "Haloalcanos",
    },
    "tox_sulfonamides": {
        "ru": "Сульфонамиды",
        "uk": "Сульфонаміди",
        "en": "Sulfonamides",
        "es": "Sulfonamidas",
    },
    "tox_phosphates": {
        "ru": "Фосфаты",
        "uk": "Фосфати",
        "en": "Phosphates",
        "es": "Fosfatos",
    },
    "tox_alert": {
        "ru": "Предупреждение",
        "uk": "Попередження",
        "en": "Alert",
        "es": "Alerta",
    },
    "tox_risk_high": {
        "ru": "ВЫСОКИЙ",
        "uk": "ВИСОКИЙ",
        "en": "HIGH",
        "es": "ALTO",
    },
    "tox_risk_medium": {
        "ru": "СРЕДНИЙ",
        "uk": "СЕРЕДНІЙ",
        "en": "MEDIUM",
        "es": "MEDIO",
    },
    "tox_risk_low": {
        "ru": "НИЗКИЙ",
        "uk": "НИЗЬКИЙ",
        "en": "LOW",
        "es": "BAJO",
    },

    # ═══════════════════════════════════════════════════════════════════════════
    # ЯМР
    # ═══════════════════════════════════════════════════════════════════════════
    "nmr_methyl": {
        "ru": "CH3 (метил)",
        "uk": "CH3 (метил)",
        "en": "CH3 (methyl)",
        "es": "CH3 (metilo)",
    },
    "nmr_methylene": {
        "ru": "CH2 (метилен)",
        "uk": "CH2 (метилен)",
        "en": "CH2 (methylene)",
        "es": "CH2 (metileno)",
    },
    "nmr_methine": {
        "ru": "CH (метин)",
        "uk": "CH (метин)",
        "en": "CH (methine)",
        "es": "CH (metino)",
    },
    "nmr_allylic": {
        "ru": "Аллильный CH",
        "uk": "Алільний CH",
        "en": "Allylic CH",
        "es": "CH alílico",
    },
    "nmr_alpha_carbonyl": {
        "ru": "CH рядом с C=O",
        "uk": "CH поруч з C=O",
        "en": "CH alpha to C=O",
        "es": "CH alfa a C=O",
    },
    "nmr_alkyne": {
        "ru": "Ацетиленовый CH",
        "uk": "Ацетиленовий CH",
        "en": "Acetylenic CH",
        "es": "CH acetilénico",
    },
    "nmr_methoxy": {
        "ru": "OCH3 (метокси)",
        "uk": "OCH3 (метокси)",
        "en": "OCH3 (methoxy)",
        "es": "OCH3 (metoxi)",
    },
    "nmr_alpha_oxygen": {
        "ru": "CH рядом с O",
        "uk": "CH поруч з O",
        "en": "CH alpha to O",
        "es": "CH alfa a O",
    },
    "nmr_alkene": {
        "ru": "Алкеновый CH",
        "uk": "Алкеновий CH",
        "en": "Alkenic CH",
        "es": "CH alquénico",
    },
    "nmr_aromatic": {
        "ru": "Ароматический CH",
        "uk": "Ароматичний CH",
        "en": "Aromatic CH",
        "es": "CH aromático",
    },
    "nmr_aldehyde": {
        "ru": "Альдегидный CH",
        "uk": "Альдегідний CH",
        "en": "Aldehyde CH",
        "es": "CH aldehído",
    },
    "nmr_carboxylic_acid": {
        "ru": "COOH (карбоновая кислота)",
        "uk": "COOH (карбонова кислота)",
        "en": "COOH (carboxylic acid)",
        "es": "COOH (ácido carboxílico)",
    },
    "nmr_phenol": {
        "ru": "OH (фенол)",
        "uk": "OH (фенол)",
        "en": "OH (phenol)",
        "es": "OH (fenol)",
    },
    "nmr_amine": {
        "ru": "NH (амин)",
        "uk": "NH (амін)",
        "en": "NH (amine)",
        "es": "NH (amina)",
    },
    "nmr_cooh": {
        "ru": "COOH (карбоновая кислота)",
        "uk": "COOH (карбонова кислота)",
        "en": "COOH (carboxylic acid)",
        "es": "COOH (ácido carboxílico)",
    },
    "nmr_ch3_near_carbonyl": {
        "ru": "CH3 соседний с C=O (ацетил/уксусная кислота)",
        "uk": "CH3 сусідній з C=O (ацетил/оцтова кислота)",
        "en": "CH3 adjacent to C=O (acetyl/acetic acid)",
        "es": "CH3 adyacente a C=O (acetil/ácido acético)",
    },
    "nmr_ch3_aliphatic": {
        "ru": "CH3 (алифатический метил)",
        "uk": "CH3 (аліфатичний метил)",
        "en": "CH3 (aliphatic methyl)",
        "es": "CH3 (metilo alifático)",
    },
    "nmr_carbonyl_ester": {
        "ru": "C=O (сложный эфир)",
        "uk": "C=O (складний естер)",
        "en": "C=O (ester)",
        "es": "C=O (éster)",
    },
    "nmr_carbonyl_amide": {
        "ru": "C=O (амид)",
        "uk": "C=O (амід)",
        "en": "C=O (amide)",
        "es": "C=O (amida)",
    },
    "nmr_carbonyl_ketone": {
        "ru": "C=O (кетон, альдегид)",
        "uk": "C=O (кетон, альдегід)",
        "en": "C=O (ketone, aldehyde)",
        "es": "C=O (cetona, aldehído)",
    },
    "nmr_nitrile": {
        "ru": "C≡N (нитрил)",
        "uk": "C≡N (нітрил)",
        "en": "C≡N (nitrile)",
        "es": "C≡N (nitrilo)",
    },
    "nmr_carbonyl_acid_ester": {
        "ru": "C=O (карбоновая кислота / сложный эфир)",
        "uk": "C=O (карбонова кислота / складний естер)",
        "en": "C=O (carboxylic acid / ester)",
        "es": "C=O (ácido carboxílico / éster)",
    },
    "nmr_carbonyl_ketone_aldehyde": {
        "ru": "C=O (кетон / альдегид)",
        "uk": "C=O (кетон / альдегід)",
        "en": "C=O (ketone / aldehyde)",
        "es": "C=O (cetona / aldehído)",
    },
    "nmr_aromatic_c": {
        "ru": "Ароматический C",
        "uk": "Ароматичний C",
        "en": "Aromatic C",
        "es": "C aromático",
    },
    "nmr_ch_methine_quaternary": {
        "ru": "CH (метин) или четвертичный углерод",
        "uk": "CH (метин) або четвертинний вуглець",
        "en": "CH (methine) or quaternary carbon",
        "es": "CH (metino) o carbono cuaternario",
    },

    # ═══════════════════════════════════════════════════════════════════════════
    # Ретросинтез
    # ═══════════════════════════════════════════════════════════════════════════
    "retro_note": {
        "ru": "(*) — точка присоединения, где был разорван связь",
        "uk": "(*) — точка приєднання, де був розірваний зв'язок",
        "en": "(*) — attachment point where bond was broken",
        "es": "(*) — punto de unión donde se rompió el enlace",
    },

    # ═══════════════════════════════════════════════════════════════════════════
    # Статусы библиотеки
    # ═══════════════════════════════════════════════════════════════════════════
    "lib_saved": {
        "ru": "сохранено",
        "uk": "збережено",
        "en": "saved",
        "es": "guardado",
    },
    "lib_updated": {
        "ru": "обновлено",
        "uk": "оновлено",
        "en": "updated",
        "es": "actualizado",
    },
    "lib_deleted": {
        "ru": "удалено",
        "uk": "видалено",
        "en": "deleted",
        "es": "eliminado",
    },

    # ═══════════════════════════════════════════════════════════════════════════
    # Индексная страница RDKit API
    # ═══════════════════════════════════════════════════════════════════════════
    "rdkit_index_title": {
        "ru": "Расширенный RDKit API",
        "uk": "Розширений RDKit API",
        "en": "Enhanced RDKit API",
        "es": "API RDKit Mejorado",
    },
    "rdkit_index_adme": {
        "ru": "ADME и Drug-likeness",
        "uk": "ADME та Drug-likeness",
        "en": "ADME & Drug-likeness",
        "es": "ADME y Similitud a Fármaco",
    },
    "rdkit_index_analysis": {
        "ru": "Анализ",
        "uk": "Аналіз",
        "en": "Analysis",
        "es": "Análisis",
    },
    "rdkit_index_reactions": {
        "ru": "Реакции и конверсия",
        "uk": "Реакції та конверсія",
        "en": "Reactions & Conversion",
        "es": "Reacciones y Conversión",
    },
    "rdkit_index_3d": {
        "ru": "3D и визуализация",
        "uk": "3D та візуалізація",
        "en": "3D & Visualization",
        "es": "3D y Visualización",
    },
    "rdkit_index_batch": {
        "ru": "Пакетная обработка и библиотека",
        "uk": "Пакетна обробка та бібліотека",
        "en": "Batch & Library",
        "es": "Lotes y Biblioteca",
    },
    "rdkit_index_comparison": {
        "ru": "Сравнение",
        "uk": "Порівняння",
        "en": "Comparison",
        "es": "Comparación",
    },
}


def get_text(key: str, lang: str = DEFAULT_LANGUAGE) -> str:
    """
    Получить перевод по ключу и коду языка.
    Если перевод не найден — возвращает ключ.
    """
    if lang not in SUPPORTED_LANGUAGES:
        lang = DEFAULT_LANGUAGE
    return TRANSLATIONS.get(key, {}).get(lang, key)


def get_lang(
    lang: str = Query(DEFAULT_LANGUAGE, description="Язык ответа: ru, uk, en, es")
) -> str:
    """
    FastAPI dependency для получения языка из query-параметра.
    """
    if lang not in SUPPORTED_LANGUAGES:
        return DEFAULT_LANGUAGE
    return lang
