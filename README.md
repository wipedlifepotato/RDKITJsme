# RadoniumAPI 

Веб-приложение для анализа молекул, хемоинформатики и планирования синтеза.
---
<img width="791" height="924" alt="image" src="https://github.com/user-attachments/assets/ed4245ff-a08c-46aa-b7d8-1de20938fd69" />
<img width="768" height="913" alt="image" src="https://github.com/user-attachments/assets/d7cb153a-9edf-4204-8830-cc1cec967a36" />
<img width="1920" height="941" alt="image" src="https://github.com/user-attachments/assets/a94cb32c-92f2-4761-a809-4d436647b91b" />

---

(periodatic table can have on commit f18b36e639055d813ae5641da31cd98d4b4bd3f0 normal biological/toxicity information, like barium is toxic, but in this version is not toxic, so recheck all information. also for on 30.09.2026 some organic things low if even it's unknown by rdkit/another places/pubchem, like unknown molecule, but it's can be harmful. Or some molecules is harmful in big concentration like acetic acid. So, all information for a now give as if but no warranty no garranty)

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
- типичные минералы

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

Для production-развёртывания используйте nginx как reverse proxy.

### Конфигурация nginx

### Пример конфигурации
nginx.conf file. Дальше уже в HTML файле прописана логика если видно, что не с локалхоста запущено

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
config.ini или аргументы командной строки. Приоритет у аргументов CLI: `-h/--host`, `-p/--port`, `-s/--socks`, `-d/--db`. Если аргумент не передан, берётся значение из config.ini.

## Лицензия
WTF
Подробнее: [LICENSE](LICENSE)

