Feature: Ensuring that Lettuce works, and the Perseq public site is loading
In order to make sure that the test system works 
As a developer
I open the perseq site using lettuce

Scenario: Opening the PerSeq site
	Given I visit the url "http://perseqer.appspot.com"
	When I look around
	I should see "PerSeq" anywhere in the page
