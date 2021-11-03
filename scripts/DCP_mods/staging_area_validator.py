"""
Runs a pre-check of a metadata files to identify issues that 
might cause the snapshot or indexing processes to fail.
Adapted from https://github.com/DataBiosphere/hca-import-validation/blob/main/hca/staging_area_validator.py
"""
import json
import os
import requests
import sys
import uuid
from typing import (
    MutableMapping,
    TypeVar,
)
from jsonschema import (
    Draft7Validator,
    FormatChecker,
    ValidationError,
)


T = TypeVar('T')
JSON = MutableMapping[str, T]
date_format = '%Y-%m-%dT%H:%M:%S.%fZ'


def validate_files(path: str, reporting) -> None:
    directory = reporting['dataset']
    validate_file_fn = globals()[f'validate_{path}_file']
    for r, d, f in os.walk(os.path.join(directory, path)):
        for file_name in f:
            full_path = os.path.join(r, file_name)
            validate_file_fn(full_path, reporting)


def validate_links_file(file_name, reporting) -> None:
    type_file = '/'.join(file_name.split('/')[1:])
    # Expected syntax: links/{bundle_uuid}_{version}_{project_uuid}.json
    if file_name.count('_') != 2:
        reporting['file_errors'].add(type_file)
        write_error("ERROR: {} has {} _'s in name, expecting 2".format(type_file, str(file_name.count('_'))))
    if not file_name.endswith('.json'):
        reporting['file_errors'].add(type_file)
        write_error('ERROR: {} does not end in .json'.format(type_file))
    _, _, project_uuid = file_name[:-5].split('_')
    file_json = json.load(open(file_name))
    validate_file_json(file_json, file_name, reporting)
    for link in file_json['links']:
        link_type = link['link_type']
        if link_type == 'process_link':
            add_metadata_file(entity_id=link['process_id'],
                                   entity_type=link['process_type'],
                                   project_uuid=project_uuid,
                                   category='process',
                                   reporting=reporting)
            for category in ('input', 'output', 'protocol'):
                for entity in link[f'{category}s']:
                    if entity.get(f'{category}_type'):
                        entity_type = entity[f'{category}_type']
                    else:
                        reporting['file_errors'].add(type_file)
                        write_error('ERROR: {}, {} missing {}_type'.format(type_file, category, category))
                        entity_type = 'missing-type'
                    if entity.get(f'{category}_id'):
                        entity_id = entity[f'{category}_id']
                    else:
                        reporting['file_errors'].add(type_file)
                        write_error('ERROR: {}, {} missing {}_id'.format(type_file, category, category))
                        entity_id = 'missing-id'
                    add_metadata_file(entity_id=entity_id,
                                           entity_type=entity_type,
                                           project_uuid=project_uuid,
                                           category=category,
                                           reporting=reporting)
        elif link_type == 'supplementary_file_link':
            if link['entity']['entity_type'] != 'project':
                reporting['file_errors'].add(type_file)
                write_error('ERROR: {}, supplementary_file not linked to entity_type project'.format(type_file))
            if link['entity']['entity_id'] != project_uuid:
                reporting['file_errors'].add(type_file)
                write_error('ERROR: {}, supplementary_file not linked to project entity_id {}'.format(type_file, project_uuid))
            for entity in link['files']:
                entity_type = entity['file_type']
                entity_id = entity['file_id']
                add_metadata_file(entity_id=entity_id,
                                   entity_type=entity_type,
                                   project_uuid=project_uuid,
                                   category='supplementary',
                                   reporting=reporting)
    if project_uuid not in reporting['metadata_files']:
        add_metadata_file(entity_id=project_uuid,
                           entity_type='project',
                           project_uuid=project_uuid,
                           category='project',
                           reporting=reporting)


def add_metadata_file(entity_id: str,
                      entity_type: str,
                      project_uuid: str,
                      category: str,
                      reporting
                      ) -> None:
    try:
        file = reporting['metadata_files'][entity_id]
    except KeyError:
        reporting['metadata_files'][entity_id] = {
            'name': None,
            'entity_id': entity_id,
            'entity_type': entity_type,
            'version': None,
            'project': {project_uuid},
            'category': {category},
            'found_metadata': False,
            'found_descriptor': False,
        }
    else:
        file['project'].add(project_uuid)
        file['category'].add(category)


def validate_metadata_file(file_name, reporting) -> None:
    type_file = '/'.join(file_name.split('/')[1:])
    # Expected syntax: metadata/{metadata_type}/{metadata_id}_{version}.json
    metadata_type, metadata_file = file_name.split('/')[-2:]
    if metadata_file.count('_') != 1:
        reporting['file_errors'].add(type_file)
        write_error("ERROR: {} has {} _'s in name, expecting 1".format(type_file, str(metadata_file.count('_'))))
    if not metadata_file.endswith('.json'):
        reporting['file_errors'].add(type_file)
        write_error('ERROR: {} does not end in .json'.format(type_file))
    metadata_id, metadata_version = metadata_file.split('_')
    file_json = json.load(open(file_name))
    validate_file_json(file_json, file_name, reporting)
    if file_json.get('provenance'):
        if metadata_id != file_json['provenance']['document_id']:
            reporting['file_errors'].add(type_file)
            write_error('ERROR: {} id in file name is {} but does not match provenance.document_id {}'.format(type_file, metadata_id, file_json['provenance']['document_id']))
    if reporting['metadata_files'].get(metadata_id):
        metadata_file = reporting['metadata_files'][metadata_id]
        metadata_file['name'] = file_name
        metadata_file['version'] = metadata_version
        metadata_file['found_metadata'] = True
        if metadata_type.endswith('_file'):
            metadata_file['data_file_name'] = file_json['file_core']['file_name']
    else:
        reporting['extra_files'].append(type_file)


def validate_descriptors_file(file_name, reporting) -> None:
    type_file = '/'.join(file_name.split('/')[1:])
    # Expected syntax: descriptors/{metadata_type}/{metadata_id}_{version}.json
    metadata_type, metadata_file = file_name.split('/')[-2:]
    if metadata_file.count('_') != 1:
        reporting['file_errors'].add(metadata_file)
        write_error("ERROR: {} has {} _'s in name, expecting 1".format(metadata_file, str(metadata_file.count('_'))))
    if not metadata_file.endswith('.json'):
        reporting['file_errors'].add(metadata_file)
        write_error('ERROR: {} does not end in .json'.format(metadata_file))
    metadata_id, metadata_version = metadata_file.split('_')
    file_json = json.load(open(file_name))
    validate_file_json(file_json, file_name, reporting)
    file_name = file_json['file_name']
    reporting['names_to_id'][file_name] = metadata_id
    if reporting['metadata_files'].get(metadata_id):
        metadata_file = reporting['metadata_files'][metadata_id]
        metadata_file['found_descriptor'] = True
        metadata_file['crc32c'] = file_json['crc32c']
        version = metadata_file['version']
        if version != metadata_version:
            reporting['file_errors'].add(metadata_file)
            write_error('ERROR: {} version in file name is {} but does not match corresponding metadata file version {}'.format(metadata_file, metadata_version, version))
    else:
        reporting['extra_files'].append(type_file)


def validate_file_json(file_json: JSON, file_name: str, reporting) -> None:
    type_file = '/'.join(file_name.split('/')[1:])
    if file_json.get('describedBy'):
        try:
            schema = download_schema(file_json['describedBy'])
        except json.decoder.JSONDecodeError as e:
            reporting['file_errors'].add(type_file)
            write_error('ERROR: {}, failed to parse schema JSON from {}'.format(type_file, file_json['describedBy']))
    else:
        reporting['file_errors'].add(type_file)
        write_error('ERROR: {}, describedBy is required'.format(type_file))
    if 'schema' in locals():
        validator = Draft7Validator(schema)
        errors = validator.iter_errors(file_json)
        for e in errors:
            reporting['file_errors'].add(type_file)
            write_error('ERROR: {}, {}'.format(type_file, e.message))


def download_schema(schema_url: str) -> JSON:
    response = requests.get(schema_url, allow_redirects=False)
    response.raise_for_status()
    return response.json()


def check_results(reporting):
    for metadata_id, metadata_file in reporting['metadata_files'].items():
        if not metadata_file['found_metadata']:
            reporting['file_errors'].add(metadata_file['entity_id'])
            write_error('ERROR: {} {} metadata file not found'.format(metadata_file['entity_type'], metadata_file['entity_id']))
        if metadata_file['entity_type'].endswith('_file') and not metadata_file['found_descriptor']:
            reporting['file_errors'].add(metadata_file['entity_id'])
            write_error('ERROR: {} {} descriptor file not found'.format(metadata_file['entity_type'], metadata_file['entity_id']))


def write_error(message):
    with open('DCP_outs/validation_errors.txt', 'a') as outfile:
        outfile.write(message + '\n')


def dcp_validation(dataset):
    write_error('validating ' + dataset)
    reporting = {
        'dataset': dataset,
        'names_to_id': {},
        'metadata_files': {},
        'file_errors': set(),
        'extra_files': [],
    }
    write_error('validating links')
    validate_files('links', reporting)
    write_error('validating metadata')
    validate_files('metadata', reporting)
    write_error('validating descriptors')
    validate_files('descriptors', reporting)
    check_results(reporting)

    for file_name in reporting['extra_files']:
        reporting['file_errors'].add(file_name)
        write_error('ERROR: {}, not a part of a subgraph in links files'.format(file_name))

    if not reporting['file_errors']:
        write_error('No DCP schema validation errors')

    return(len(reporting['file_errors']))
