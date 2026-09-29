# RadoniumAPI 

Веб-приложение для анализа молекул, хемоинформатики и планирования синтеза.
---
<img width="791" height="924" alt="image" src="https://github.com/user-attachments/assets/ed4245ff-a08c-46aa-b7d8-1de20938fd69" />
<img width="768" height="913" alt="image" src="https://github.com/user-attachments/assets/d7cb153a-9edf-4204-8830-cc1cec967a36" />
---
(periodatic table can have on commit f18b36e639055d813ae5641da31cd98d4b4bd3f0 normal biological/toxicity information, like barium is toxic, but in this version is not toxic, so recheck all information)

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

## Быстрый старт

```bash
# Клонирование репозитория
git clone https://github.com/yourusername/chem.git
cd chem

# Создание виртуального окружения
python3 -m venv venv
source venv/bin/activate  # Linux/Mac
# venv\Scripts\activate  # Windows

# Установка зависимостей
pip install -r reqs.txt

# Запуск приложения
python3 main.py
```

Откройте браузер: http://localhost:1235

## Nginx

Для production-развёртывания используйте nginx как reverse proxy.

### Конфигурация nginx

Файл `nginx.conf` содержит пример конфигурации. Перед использованием необходимо изменить:

1. **Пути к файлам:**
   - `/opt/chem/RDKITJsme/examples_clients/` — путь к клиентским файлам
   - `/opt/chem/RDKITJsme/StaticFiles/` — путь к статическим файлам
   - `/var/www/html` — путь к статическим файлам по умолчанию

2. **Порт приложения:**
   - В блоке `location /rdkit/` указан `proxy_pass http://127.0.0.1:1235/;`
   - Измените порт `1235` на порт, который использует ваше приложение (по умолчанию `1235`)

3. **SSL (опционально):**
   - Раскомментируйте строки с `ssl_certificate` и `ssl_certificate_key`
   - Укажите пути к вашим SSL-сертификатам

4. **Server name (опционально):**
   - Раскомментируйте `server_name` и укажите ваш домен

### Пример конфигурации
nginx.conf file

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

## Использование

### Веб-интерфейс

1. Ввод SMILES: Введите SMILES-строку в поле ввода
2. Редактирование: Используйте JSME-редактор для рисования молекул
3. Анализ: Переключайтесь между вкладками (Basic, ADME, Analysis, Advanced, 3D, NMR, Toxicity, Retrosynthesis, Conformers, Export)
4. Справка: Нажмите кнопку "Help" для объяснения всех дескрипторов

### Примеры SMILES

| Молекула | SMILES |
|----------|--------|
| Аспирин | `CC(=O)Oc1ccccc1C(=O)O` |
| Кофеин | `CN1C=NC2=C1C(=O)N(C(=O)N2C)C` |
| Этанол | `CCO` |
| Бензол | `c1ccccc1` |
| Глюкоза | `C([C@@H]1[C@H]([C@@H]([C@H](C(O1)O)O)O)O)O` |
| Кетозы длиной цепи глюкозы (для проверки стереоизомерии, функционала) (фруктоза, сорбоза, псикоза, тагатоза)  | `O=C(CO)C(O)C(O)C(O)CO` |

## API Endpoints

### Основные

| Метод | Путь | Описание |
|-------|------|----------|
| GET | `/` | Главная страница (JSME-редактор) |
| GET | `/api/get_properties` | Базовые свойства молекулы |
| GET | `/api/convert` | Конвертация форматов |
| GET | `/api/similarity` | Сравнение двух молекул |
| GET | `/api/get_3d_sdf` | 3D-структура в SDF |
| GET | `/api/get_name` | Название молекулы (IUPAC) |
| GET | `/api/get_chiral` | Хиральные центры |
| GET | `/api/substructure_search` | Поиск подструктур |
| GET | `/api/render` | Рендеринг 2D |

### Расширенные (RDKit)

| Метод | Путь | Описание |
|-------|------|----------|
| GET | `/rdkit/api/adme` | ADME-свойства, QED, SA Score |
| GET | `/rdkit/api/functional_groups` | Поиск функциональных групп |
| GET | `/rdkit/api/validate` | Валидация SMILES |
| GET | `/rdkit/api/pka` | Оценка pKa |
| GET | `/rdkit/api/descriptors` | Все дескрипторы (200+) |
| GET | `/rdkit/api/reaction` | Применение реакций |
| GET | `/rdkit/api/by_inchikey` | Поиск по InChIKey |
| GET | `/rdkit/api/3d` | 3D-структура |
| GET | `/rdkit/api/nmr` | Оценка ЯМР (1H и 13C) |
| GET | `/rdkit/api/toxicity` | Оценка токсичности |
| GET | `/rdkit/api/retrosynthesis` | Ретросинтетический анализ |
| GET | `/rdkit/api/conformers` | Генерация конформеров |
| GET | `/rdkit/api/export_svg` | Экспорт в SVG |
| GET | `/rdkit/api/ms` | Масс-спектрометрия |
| POST | `/rdkit/api/batch` | Пакетный анализ |
| GET | `/rdkit/api/library` | Библиотека молекул |
| POST | `/rdkit/api/library/save` | Сохранение в библиотеку |
| DELETE | `/rdkit/api/library/{smiles}` | Удаление из библиотеки |
| POST | `/rdkit/api/compare` | Сравнение молекул |
| GET | `/rdkit/api/render_highlighted` | Рендеринг с подсветкой |

### Документация API

- Swagger UI: http://localhost:1235/docs
- ReDoc: http://localhost:1235/redoc
- OpenAPI JSON: http://localhost:1235/openapi.json

## Конфигурация

### Прокси (SOCKS5)

```ini
# config.ini
[proxy]
address = 127.0.0.1:9050
```

### Порт и хост
config.ini или аргументы командной строки. Приоритет у аргументов CLI: `-h/--host`, `-p/--port`, `-s/--socks`, `-d/--db`. Если аргумент не передан, берётся значение из config.ini.

## Лицензия

Personal & Educational Use License

- Бесплатно: личное использование, образование, некоммерческие исследования, Open Source
- Платно: коммерческое использование, SaaS, встраивание в платные продукты

Для получения коммерческой лицензии свяжитесь: issues github

Подробнее: [LICENSE](LICENSE)

