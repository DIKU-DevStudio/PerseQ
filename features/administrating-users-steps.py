from lettuce import step

@step(u'Given I am on the user management page')
def given_i_am_on_the_user_management_page(step):
    assert False, 'This step must be implemented'

@step(u'Then I should see a list of users')
def then_i_should_see_a_list_of_users(step):
    assert False, 'This step must be implemented'

@step(u'And there is a user called "([^"]*)" in the list')
def and_there_is_a_user_called_group1_in_the_list(step, group1):
    assert False, 'This step must be implemented'

@step(u'When I click the delete button for "([^"]*)"')
def when_i_click_the_delete_button_for_group1(step, group1):
    assert False, 'This step must be implemented'

@step(u'Then there should be no user called "([^"]*)" in the list')
def then_there_should_be_no_user_called_group1_in_the_list(step, group1):
    assert False, 'This step must be implemented'

@step(u'Given I am logged in as a normal user')
def given_i_am_logged_in_as_a_normal_user(step):
    assert False, 'This step must be implemented'

@step(u'And I am on the user management page')
def and_i_am_on_the_user_management_page(step):
    assert False, 'This step must be implemented'

@step(u'Then I should see "([^"]*)" anywhere on the page')
def then_i_should_see_group1_anywhere_on_the_page(step, group1):
    assert False, 'This step must be implemented'