import ee
cloud_project = 'ee-tolbico2024'

try:
        ee.Initialize(project=cloud_project)
except:
        ee.Authenticate()
        ee.Initialize(project=cloud_project)
import ee
cloud_project = 'ee-tolbico2024'

try:
        ee.Initialize(project=cloud_project)
except:
        ee.Authenticate()
        ee.Initialize(project=cloud_project)
