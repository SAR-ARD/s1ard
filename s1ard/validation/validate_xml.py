from lxml import etree
import sys
import os


def validate(xml_path, xsd_path):
    message = ''
    try:
        xsd = etree.parse(xsd_path)
        xml = etree.parse(xml_path)
        schema = etree.XMLSchema(xsd)
        schema.validate(xml)
        message += str(schema.error_log)
    except Exception as e:
        message += str(e)
    finally:
        if message != '':
            print('{} is NOT VALID:\n{}'.format(xml_path, message))
            return 1
        else:
            print('{} is VALID!'.format(xml_path))
            return 0


if __name__ == '__main__':
    xml_path = sys.argv[1]
    xsd_path = sys.argv[2]
    if not xml_path.endswith('.xml') or not xsd_path.endswith('.xsd'):
        raise RuntimeError('Invalid input! Correct usage of this script:\n\n'
                           'python validate_xml.py your_xml_file.xml your_xsd_file.xsd\n')
    if any([not os.path.isfile(x) for x in [xml_path, xsd_path]]):
        raise FileNotFoundError('Input file(s) not found!')
    
    ret = validate(xml_path=xml_path, xsd_path=xsd_path)
    sys.exit(ret)
