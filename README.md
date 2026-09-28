# SMILES API RDKit

Веб-приложение для анализа молекул, хемоинформатики и планирования синтеза.

## Возможности

### Медицинская химия
- Базовые свойства: Молекулярная масса, LogP, TPSA, HBD, HBA
- Лекарственность: QED, Правило Липинского, SA Score
- Функциональные группы: Поиск 25+ типов функциональных групп

### Анализ и реакции
- pKa оценка: Кислотные и основные центры
- 200+ дескрипторов: Полный набор дескрипторов RDKit
- Реакции: Применение SMIRKS-реакций к молекулам

### ЯМР и спектроскопия
- 1H ЯМР: Эвристическая оценка химических сдвигов протонов
- 13C ЯМР: Эвристическая оценка химических сдвигов углерода
- Масс-спектрометрия: Изотопное распределение

### Токсичность (ADMET)
- RDKit alerts

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

Откройте браузер: http://localhost:8000

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
   - Измените порт `1235` на порт, который использует ваше приложение (по умолчанию `8000`)

3. **SSL (опционально):**
   - Раскомментируйте строки с `ssl_certificate` и `ssl_certificate_key`
   - Укажите пути к вашим SSL-сертификатам

4. **Server name (опционально):**
   - Раскомментируйте `server_name` и укажите ваш домен

### Пример конфигурации

```nginx
server {
    listen 80;
    server_name your-domain.com;

    # Все запросы к /rdkit/ (UI + эндпоинты API)
    location /rdkit/ {
        proxy_pass http://127.0.0.1:8000/;
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;
    }

    # Документация FastAPI
    location /docs { proxy_pass http://127.0.0.1:8000; }
    location /redoc { proxy_pass http://127.0.0.1:8000; }
    location /openapi.json { proxy_pass http://127.0.0.1:8000; }

    # Статические файлы
    location /StaticFiles/ {
        alias /path/to/your/StaticFiles/;
        expires 7d;
        add_header Cache-Control "public, no-transform";
    }
}
```

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

- Swagger UI: http://localhost:8000/docs
- ReDoc: http://localhost:8000/redoc
- OpenAPI JSON: http://localhost:8000/openapi.json

## Конфигурация

### Прокси (SOCKS5)

```ini
# config.ini
[proxy]
address = 127.0.0.1:9050
```

### Порт и хост

```python
# main.py
uvicorn.run(app, host="0.0.0.0", port=8000)
```

## Лицензия

Personal & Educational Use License

- Бесплатно: личное использование, образование, некоммерческие исследования, Open Source
- Платно: коммерческое использование, SaaS, встраивание в платные продукты

Для получения коммерческой лицензии свяжитесь: issues github

Подробнее: [LICENSE](LICENSE)

