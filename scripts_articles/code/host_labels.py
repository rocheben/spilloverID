"""Display labels for candidate hosts, shared by the figure scripts.

Common name first, scientific name in italics in parentheses, so that the
register is consistent across panels (cattle / wild ungulates / pigs are
common terms, Pteropus and Peromyscus are taxa and are italicised).
"""

DISPLAY = {
    # Sin Nombre
    'P. maniculatus': r'Deer mouse ($\it{P.\ maniculatus}$)',
    'P. truei':       'Pi\u00f1on mouse ($\\it{P.\\ truei}$)',
    'Microtus':       r'Voles ($\it{Microtus}$)',
    'Other rodents':  'Other rodents',
    # Nipah
    'Pteropus':       r'Flying foxes ($\it{Pteropus}$)',
    'Pigs':           'Pigs',
    'Equines':        'Equines',
    'Dogs':           'Dogs',
    # Rift Valley fever
    'Wild ungulates':  'Wild ungulates',
    'Cattle':          'Cattle',
    'Camels':          'Camels',
    'Small ruminants': 'Small ruminants',
}


def label(host):
    return DISPLAY.get(host, host)


def labels(hosts):
    return [label(h) for h in hosts]
