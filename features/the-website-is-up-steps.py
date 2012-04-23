from salad.steps.everything import *

@step('I should see "([^"]*)" somewhere')
def i_should_see_group1_somewhere(step, group1):
	assert False, 'This step must be implemented'