<div align="center">

# RadoniumAPI

<p>
  <a href="#english"><b>English</b></a> &nbsp;&nbsp;|&nbsp;&nbsp;
  <a href="#russian"><b>Русский</b></a>
</p>

</div>

---

<div id="english"></div>

# RadoniumAPI

A web application for molecular analysis, chemoinformatics, and synthesis planning.
---
<img width="791" height="924" alt="image" src="https://github.com/user-attachments/assets/ed4245ff-a08c-46aa-b7d8-1de20938fd69" />
<img width="768" height="913" alt="image" src="https://github.com/user-attachments/assets/d7cb153a-9edf-4204-8830-cc1cec967a36" />
<img width="1920" height="941" alt="image" src="https://github.com/user-attachments/assets/a94cb32c-92f2-4761-a809-4d436647b91b" />

---

*(Note: Periodic table information around commit `f18b36e639055d813ae5641da31cd98d4b4bd3f0` reflects normal biological/toxicity data—e.g., barium being treated as non-toxic in that specific version, so please recheck critical data. Also, as of September 30, 2026, some obscure or novel organic molecules might be classified as low-risk or unknown by RDKit/PubChem even if potentially hazardous, and substances like acetic acid are harmful at high concentrations. All information is provided "as is" with no warranty or guarantee.)*

## Features

### Medicinal Chemistry
- Basic Properties: Molecular weight, LogP, TPSA, HBD, HBA
- Drug-likeness: QED, Lipinski's Rule of Five, SA Score
- Functional Groups: Detection of 25+ functional group types

### Analysis & Reactions
- pKa Estimation: Acidic and basic centers
- 200+ Descriptors: Full set of RDKit descriptors

### NMR & Spectroscopy
- 1H NMR: Heuristic proton chemical shift estimation
- 13C NMR: Heuristic carbon chemical shift estimation
- Mass Spectrometry: Isotopic distribution

### Toxicity (ADMET)
- RDKit alerts
- PubChem GHS Classification

### Retrosynthesis
- Bond cleavage analysis (RDKit RECAP)

### Visualization
- 2D Rendering: With chiral center highlighting
- 3D Viewer: Interactive 3D viewer (3Dmol.js)
- Conformers: Generation of 5-10 low-energy conformers
- Paint Toolbar: Drawing over 2D renders

### Export
- SVG: Vector graphics for research publications
- SDF: 3D structures
- PNG: With annotations

### Tools
- SMILES Conversion: To InChI, InChIKey, canonical SMILES
- Substructure Search: SMARTS patterns
- Molecule Comparison: Tanimoto similarity
- Library: Save and manage molecules

### Periodic Table
- Interactive table with traditional grouping (1A, 2A, 3B... 8A)
- Color coding: alkali, alkaline earth, transition, post-transition, metalloids, nonmetals, halogens, noble gases, lanthanides, actinides
- Element click: spdf electron configuration, valence electron spins, ion color, toxicity, bioavailability
- Dobereiner-Mendeleev triads rule
- Category filters
- Multilingual: Russian, English, Ukrainian, Spanish
- Access: http://localhost:1235/periodatic_table
- Typical minerals

## Quick Start

```bash
# Clone the repository
git clone https://github.com/wipedlifepotato/RadoniumAPI 
cd RadoniumAPI

# Create a virtual environment
python3 -m venv venv
source venv/bin/activate  # Linux/Mac
# venv\Scripts\activate  # Windows

# Install dependencies
python3 -m pip install -r reqs.txt

# Run the application
python3 main.py
```

Open your browser: http://localhost:1235

## Nginx

For production deployment, use Nginx as a reverse proxy.

### Nginx Configuration
`nginx.conf` file. Additional routing logic is handled inside the HTML if launched outside of localhost.

### Systemd Service

Create the file `/etc/systemd/system/chem.service`:

```ini
[Unit]
Description=Chem SMILES API
After=network.target

[Service]
Type=simple
User=your-user
WorkingDirectory=/path/to/chem
ExecStart=/path/to/chem/venv/bin/python3 main.py
Restart=always

[Install]
WantedBy=multi-user.target
```

Commands to run:
```bash
sudo systemctl enable chem
sudo systemctl start chem
sudo systemctl status chem
```

## Documentation

- Swagger UI: http://localhost:1235/docs
- ReDoc: http://localhost:1235/redoc

## Configuration

### Port and Host
Configured via `config.ini` or command-line arguments. CLI arguments take priority: `-h/--host`, `-p/--port`, `-s/--socks`, `-d/--db`. If an argument is omitted, the value from `config.ini` is used.

## License
[LICENSE](LICENSE)

---

<div id="russian"></div>

# RadoniumAPI (Русская версия)

Веб-приложение для анализа молекул, хемоинформатики и планирования синтеза.
---
<img width="791" height="924" alt="image" src="https://github.com/user-attachments/assets/ed4245ff-a08c-46aa-b7d8-1de20938fd69" />
<img width="768" height="913" alt="image" src="https://github.com/user-attachments/assets/d7cb153a-9edf-4204-8830-cc1cec967a36" />
<img width="1920" height="941" alt="image" src="https://github.com/user-attachments/assets/a94cb32c-92f2-4761-a809-4d436647b91b" />

---

*(Примечание: периодическая таблица на коммите `f18b36e639055d813ae5641da31cd98d4b4bd3f0` может содержать упрощённые данные по токсичности, например, барий в ней помечен как нетоксичный, поэтому перепроверяйте информацию. Кроме того, по состоянию на 30.09.2026 некоторые органические соединения могут распознаваться RDKit / PubChem как неизвестные или малоопасные, хотя потенциально токсичны, либо опасны в больших концентрациях, как уксусная кислота. Вся информация предоставляется «как есть» без каких-либо гарантий.)*

## Возможности

### Медицинская химия
- Базовые свойства: Молекулярная масса, LogP, TPSA, HBD, HBA
- Лекарственность: QED, Правило Липинского, SA Score
- Функциональные группы: Поиск 25+ типов функциональных групп

### Анализ и реакции
- pKa оценка: Кислотные и основные центры
- 200+ дескрипторов: Полный набор дескрипторов RDKit

### ЯМР и спектроскопия
- 1H ЯМР: Эвристическая оценка химических сдвигов протонов
- 13C ЯМР: Эвристическая оценка химических сдвигов углерода
- Масс-спектрометрия: Изотопное распределение

### Токсичность (ADMET)
- RDKit alerts
- PubChem GHS Classification

### Ретросинтез
- Анализ разрывов связей (RDKit RECAP)

### Визуализация
- 2D рендеринг: С подсветкой хиральных центров
- 3D просмотр: Интерактивный 3D-вьюер (3Dmol.js)
- Конформеры: Генерация 5-10 энергетически выгодных конформеров
- Paint Toolbar: Рисование поверх 2D-рендера

### Экспорт
- SVG: Векторная графика для научных публикаций
- SDF: 3D-структуры
- PNG: С аннотациями

### Инструменты
- SMILES конвертация: В InChI, InChIKey, канонический SMILES
- Поиск подструктур: SMARTS-паттерны
- Сравнение молекул: Танимото сходство
- Библиотека: Сохранение и управление молекулами

### Периодическая таблица
- Интерактивная таблица с русской группировкой (1А, 2А, 3Б... 8А)
- Цветовая подсветка: щелочные, щёлочноземельные, переходные, постпереходные, полуметаллы, неметаллы, галогены, благородные газы, лантаноиды, актиноиды
- Клик по элементу: электронная конфигурация (spdf), спины валентных электронов, цвет иона, токсичность, биодоступность
- Правило триад (Дёберейнер — Менделеев)
- Фильтры по категориям элементов
- Мультиязычность: русский, английский, украинский, испанский
- Доступ: http://localhost:1235/periodatic_table
- Типичные минералы

## Быстрый старт

```bash
# Клонирование репозитория
git clone https://github.com/wipedlifepotato/RadoniumAPI 
cd RadoniumAPI

# Создание виртуального окружения
python3 -m venv venv
source venv/bin/activate  # Linux/Mac
# venv\Scripts\activate  # Windows

# Установка зависимостей
python3 -m pip install -r reqs.txt

# Запуск приложения
python3 main.py
```

Откройте браузер: http://localhost:1235

## Nginx

Для production-развёртывания используйте Nginx как reverse proxy.

### Конфигурация Nginx
Файл `nginx.conf`. Дальше логика прописана прямо в HTML-файле, если приложение запущено не на localhost.

### Запуск через systemd

Создайте файл `/etc/systemd/system/chem.service`:

```ini
[Unit]
Description=Chem SMILES API
After=network.target

[Service]
Type=simple
User=your-user
WorkingDirectory=/path/to/chem
ExecStart=/path/to/chem/venv/bin/python3 main.py
Restart=always

[Install]
WantedBy=multi-user.target
```

Запуск:
```bash
sudo systemctl enable chem
sudo systemctl start chem
sudo systemctl status chem
```

## Доки

- Swagger UI: http://localhost:1235/docs
- ReDoc: http://localhost:1235/redoc

## Конфигурация

### Порт и хост
`config.ini` или аргументы командной строки. Приоритет у аргументов CLI: `-h/--host`, `-p/--port`, `-s/--socks`, `-d/--db`. Если аргумент не передан, берётся значение из `config.ini`.

## Лицензия
[LICENSE](LICENSE)
