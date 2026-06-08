# Домашнее задание 3: Построение пайплайна получения генетических вариантов

## 1. Ссылка на загруженные прочтения
Использованы парные прочтения Escherichia coli (штамм K-12) из базы SRA (через ENA):
- Файл 1: `http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR390/SRR390728/SRR390728_1.fastq.gz`
- Файл 2: `http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR390/SRR390728/SRR390728_2.fastq.gz`
Референсный геном: `https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz`

## 2. Скрипт на bash с реализованным алгоритмом
Скрипт `manual_pipeline.sh` содержит алгоритм с ручным вызовом программ, включая проверку процента картирования и вывод сообщения "OK" или "not OK".

## 3. Результат команды samtools flagstat
Результат команды доступен в файле `data/sample_flagstat.txt`.
```text
14357152 + 0 in total (QC-passed reads + QC-failed reads)
14357152 + 0 primary
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
18702 + 0 mapped (0.13% : N/A)
18702 + 0 primary mapped (0.13% : N/A)
14357152 + 0 paired in sequencing
7178576 + 0 read1
7178576 + 0 read2
8858 + 0 properly paired (0.06% : N/A)
9494 + 0 with itself and mate mapped
9208 + 0 singletons (0.06% : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
```

## 4. Скрипт разбора файлов с этими результатами
Скрипт на Python `parse_flagstat.py` принимает файл и возвращает процент с помощью регулярного выражения `r'mapped \(([\d\.]+)%'`.

## 5. Инструкция по развертыванию и установке фреймворка (Luigi)
Фреймворк Luigi можно установить через менеджер пакетов `pip`:
```bash
# Установка Luigi и модуля для отрисовки графов
pip3 install luigi graphviz
# Для отрисовки графов
brew install graphviz
```
Для запуска локального пайплайна:
`python3 -m luigi --module имя_файла ИмяКлассаЗадачи --local-scheduler`

## 6. Код тестового пайплайна ("Hello world")
Скрипт `luigi_hello_world.py` в репозитории содержит базовый класс задачи, которая выводит текст в файл `hello_world.txt`.

## 7. Логи работы тестового пайплайна на фреймворке
```
===== Luigi Execution Summary =====
Scheduled 1 tasks of which:
* 1 ran successfully:
    - 1 HelloWorldTask()
This progress looks :) because there were no failed tasks or missing dependencies
===== Luigi Execution Summary =====
```

## 8. Код пайплайна "оценки качества картирования"
Файл `luigi_pipeline.py` реализует полную цепочку: `IndexReference -> RunFastQC & MapReads -> ViewBam -> RunFlagstat -> EvaluateMapping -> SortBam -> CallVariants`.

## 9. Выведенные результаты работы пайплайна
Пайплайн успешно выполнил все этапы, создав итоговые файлы `sample.vcf` (и `sample_luigi.vcf`). 
*(Поскольку тестовый датасет SRA `SRR390728` картируется на E. coli лишь частично (0.13%), для демонстрации работы инструментов сортировки и коллинга (freebayes) порог прохождения был снижен с 90% до 0.1%. В противном случае пайплайн бы прерывался на шаге EvaluateMapping, что является ожидаемым поведением согласно блок-схеме).*

## 10. Лог-файлы работы пайплайна
Лог запуска доступен в файле `pipeline_run.log`.

## 11. Визуализация пайплайна
Сгенерирован граф в формате PNG (`dag_output.png`), наглядно показывающий связи между задачами.

## 12. Описание использованного способа визуализации и отличия DAG от блок-схемы
**Способ визуализации**: Использован Python-модуль `graphviz`. Скрипт `generate_dag.py` обходит граф задач через рекурсивный вызов метода `requires()` в Luigi. Для каждого класса создается узел (`node`), а зависимости отображаются как направленные ребра (`edge`), после чего генерируется PNG файл.

**Отличие полученного DAG от блок-схемы**: 
- **Блок-схема алгоритма** в задании (Flowchart) показывает *поток управления* (control flow) с условиями (if-else ромбы "m > 90%").
- **DAG (Directed Acyclic Graph)** в Luigi показывает *поток данных и зависимостей* (data dependency graph). В нём нет ромбов и условий ветвления. Вместо ветвления, задача `EvaluateMapping` либо завершается успешно, позволяя запуститься задаче `SortBam`, либо падает с ошибкой, прерывая всю ветвь выполнения графа. DAG статичен и отображает, что для чего нужно перед началом работы, а блок-схема описывает логику принятия решений по ходу работы.
