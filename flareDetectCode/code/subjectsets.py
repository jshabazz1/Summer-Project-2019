from panoptes_client import Panoptes, Project, SubjectSet, Subject

Panoptes.connect(username='Jshabazz', password='JTzeez6062291$')


summer_project = Project.find(9645)
print(summer_project.display_name)

subject_set = SubjectSet()

# subject_set.links.project = summer_project
# subject_set.display_name = 'Flare subject set'

# subject_set.save()

project.reload()
print(project.links.subject_sets)

subject_metadata = {
    '/Users/jshabazz/Work/lightcurves/': {
        'subject_reference': 1,
        'date': '2017-01-01',
    },
    '/Users/jshabazz/Work/lightcurves/': {
        'subject_reference': 2,
        'date': '2017-01-02',
    },
    '/Users/jshabazz/Work/lightcurves/': {
        'subject_reference': 3,
        'date': '2017-01-03',
    },
}

new_subjects = []

for filename, metadata in subject_metadata.items():
    subject = Subject()

    subject.links.project = summer_project
    subject.add_location(filename)

    subject.metadata.update(metadata)

    subject.save()
    new_subjects.append(subject)