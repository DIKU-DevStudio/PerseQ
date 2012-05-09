Feature: Perseq administrators should be able to manage users
In order to keep track of user accounts which are used to track favourites, comments etc.
As an administrator
I want to be able to manage user accounts

Scenario: See a list of users
    Given I am on the user management page
    When I look around
    Then I should see a list of users

Scenario: Deleting a user
    Given I am on the user management page
    And there is a user called "deleteme" in the list
    When I click the delete button for "deleteme"
    Then there should be no user called "deleteme" in the list

Scenario: Change a user's group
    Given I am on the user management page
    And there is a user called "changemygroup" in the list

Scenario: Nothing visible if user is not an administrator
    Given I am logged in as a normal user
    And I am on the user management page
    Then I should see "You cannot administer users" anywhere on the page