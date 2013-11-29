# Copyright (C) 2011-2012 Johan Hake
#
# This file is part of Gotran.
#
# Gotran is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Gotran is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with Gotran. If not, see <http://www.gnu.org/licenses/>.

import sys, os, re
import urllib
from xml.etree import ElementTree

from collections import OrderedDict, deque, defaultdict
from gotran.common import warning, error, check_arg

from modelparameters.codegeneration import _all_keywords
from modelparameters.parameterdict import *

from mathml import MathMLBaseParser

__all__ = ["cellml2ode", "CellMLParser"]

si_unit_map = {"ampere":"A", "becquerel":"Bq", "candela":"cd", "celsius":"gradC",
               "coulomb":"C","dimensionless":"1", "farad":"F", "gram":"g",
               "gray":"Gy", "henry":"H", "hertz":"Hz", "joule":"J", "katal":"kat",
               "kelvin":"K", "kilogram":"kg", "liter":"l", "litre":"l",
               "lumen":"lm", "lux":"lx", "meter":"m", "metre":"m", "mole":"mole",
               "newton":"N", "ohm":"Omega", "pascal":"Pa", "radian":"rad",
               "second":"s", "siemens":"S", "sievert":"Sv", "steradian":"sr",
               "tesla":"T", "volt":"V", "watt":"W", "weber":"Wb", "dimensionless":"1"}

prefix_map = {"deca":"da", "hecto":"h", "kilo":"k", "mega":"M", "giga":"G",
              "tera":"T", "peta":"P", "exa":"E", "zetta":"Z", "yotta":"Y",
              "deci":"d", "centi":"c", "milli":"m", "micro":"u", "nano":"n",
              "pico":"p", "femto":"f", "atto":"a", "zepto":"z", "yocto":"y",
              None:""}

ui = "UNINITIALIZED"

class Equation(object):
    """
    Class for holding information about an Equation
    """
    def __init__(self, name, expr, used_variables):
        self.name = name
        self.expr = expr
        self.used_variables = used_variables
        self.dependent_equations = []
        self.component = None

    def check_dependencies(self, equation):
        """
        Check Equation dependencies
        """
        assert(isinstance(equation, Equation))

        if equation.name in self.used_variables:
            self.dependent_equations.append(equation)
            
    def __str__(self):
        return self.name
    
    def __repr__(self):
        return "Equation({0} = {1})".format(self.name, "".join(self.expr))

    def __eq__(self, other):
        if not isinstance(other, Equation):
            return False
        return other.name == self.name and other.component == self.component

    def __hash__(self):
        return hash(self.name+self.component.name)

class Component(object):
    """
    Class for holding information about a CellML Component
    """
    def __init__(self, name, variables, equations, state_variables=None):
        self.name = name

        self.variable_types = {}
        
        self.state_variables = OrderedDict((state, variables.pop(state, None))\
                                           for state in state_variables)
        self.variable_types.update((state, "state_variable") \
                                   for state in self.state_variables)

        self.parameters = OrderedDict((name, info) for name, info in \
                                      variables.items() if info["init"] is not None)

        self.variable_types.update((param, "parameters") \
                                   for param in self.parameters)

        self.derivatives = state_variables

        self.parent = None
        self.children = []
        
        self.dependencies = OrderedDict() 
        self.used_in = OrderedDict() 

        # Store equations
        self.sort_and_store_equations(equations)

        # Extract info
        dummy = dict(init=None, unit="1", private=True)
        self.equation_info = OrderedDict((eq.name, variables.get(eq.name, dummy))\
                                         for eq in self.equations)
        
        self.variable_types.update((eq.name, "equation") \
                                   for eq in self.equations)
        
        # Get used variables
        self.used_variables = set()
        for equation in self.equations:
            self.used_variables.update(equation.used_variables)

        # Remove dependencies on names defined by component
        self.used_variables.difference_update(\
            equation.name for equation in self.equations)

        self.used_variables.difference_update(\
            name for name in self.parameters)

        self.used_variables.difference_update(\
            name for name in self.state_variables)

    def sort_and_store_equations(self, equations):
        
        # Check internal dependencies
        for eq0 in equations:
            
            eq0.dependent_equations = []
            
            # Store component
            eq0.component = self
            for eq1 in equations:
                if eq0 == eq1:
                    continue
                eq0.check_dependencies(eq1)

        sorted_equations = []
        while equations:
            equation = equations.pop(0)
            if any(dep in equations for dep in equation.dependent_equations):
                equations.append(equation)
            else:
                sorted_equations.append(equation)

        # Store the sorted equations
        self.equations = sorted_equations


    def __hash__(self):
        return hash(self.name)

    def __str__(self):
        return self.name + "<{0}>".format(len(self.state_variables))

    def __repr__(self):
        return "Component<{0}, {1}>".format(self.name, \
                                            len(self.state_variables))

    def __eq__(self, other):
        if not isinstance(other, Component):
            return False
        
        return other.name == self.name

    def check_dependencies(self, component):
        """
        Check components dependencies
        """
        assert(isinstance(component, Component))

        if any(equation.name in self.used_variables \
               for equation in component.equations):
            dep_equations = [equation for equation in component.equations \
                             if equation.name in self.used_variables]

            # Register mutual dependencies
            self.dependencies[component] = dep_equations
            component.used_in[self] = dep_equations
            
            # Add logics for dependencies of all Equations.
            for other_equation in component.equations:
                for equation in self.equations:
                    if other_equation.name in equation.used_variables:
                        equation.dependent_equations.append(other_equation)

    def change_parameter_name(self, oldname, newname=None):
        """
        Change the name of a parameter
        Assume the name is only used locally within this component
        """
        assert(oldname in self.parameters)

        # If no newname is give we pick one based on the component name
        if newname is None:
            newname = oldname + "_" + self.name.split("_")[0]
        
        warning("Locally change parameter name '{0}' to '{1}' in "\
                "component '{2}'.".format(oldname, newname, self.name))
        
        # Update parameters
        self.parameters = OrderedDict((newname if name == oldname else \
                                       name, value) for name, value in \
                                      self.parameters.items())

        # Update equations
        for eqn in self.equations:
            while oldname in eqn.expr:
                eqn.expr[eqn.expr.index(oldname)] = newname

        return newname

    def change_state_name(self, oldname, newname=None):
        """
        Change the name of a state
        Assume the name is only used locally within this component
        """
        assert(oldname in self.state_variables)
        
        # If no newname is give we pick one based on the component name
        if newname is None:
            newname = oldname + "_" + self.name.split("_")[0]
        
        warning("Locally change state name '{0}' to '{1}' in component "\
                "'{2}'.".format(oldname, newname, self.name))
        
        # Update parameters
        self.state_variables = OrderedDict((newname if name == oldname \
                            else name, value) for name, value in \
                            self.state_variables.items())

        oldder = self.derivatives[oldname]
        newder = oldder.replace(oldname, newname)
        self.derivatives = OrderedDict((newname if name == oldname else name, \
                                        newder if value == oldder else value) \
                                       for name, value in \
                                       self.derivatives.items())
        # Update equations
        for eqn in self.equations:
            while oldname in eqn.expr:
                eqn.expr[eqn.expr.index(oldname)] = newname
            while oldder in eqn.expr:
                eqn.expr[eqn.expr.index(oldder)] = newder
            if oldder == eqn.name:
                eqn.name = newder

        old_eq_name = "d{0}_dt".format(oldname)
        new_eq_name = "d{0}_dt".format(newname)
        self.variable_types[new_eq_name]  = self.variable_types.pop(\
            old_eq_name, "equation")
        self.equation_info[new_eq_name] = self.equation_info.pop(\
            old_eq_name, dict(init=None, unit="1", private=True))

        return newname

    def change_equation_name(self, oldname, newname=None):
        
        assert(self.variable_types.get(oldname)=="equation")
        
        # If no newname is give we pick one based on the component name
        if newname is None:
            newname = oldname + "_" + self.name.split("_")[0]
        
        warning("Locally change equation name '{0}' to '{1}' in "\
                "component '{2}'.".format(oldname, newname, self.name))
        
        # Update equations
        new_equations = []
        for eqn in self.equations:
            while oldname in eqn.expr:
                eqn.expr[eqn.expr.index(oldname)] = newname

            if eqn.name == oldname:
                eqn.name = newname

            new_equations.append(eqn)

        self.equations = new_equations
        eq_info = self.equation_info.pop(oldname)
        self.equation_info[newname] = eq_info
        
        return newname

    def change_variable_name(self, oldname, newname=None):
        vartype = self.variable_types.get(oldname)

        if vartype is None:
            error("Cannot change variable name. {0} is not a variable "\
                  "in component {1}".format(oldname, self.name))
            return 
            
        if vartype == "state_variable":
            return self.change_state_name(oldname, newname)
        if vartype == "parameter":
            return self.change_parameter_name(oldname, newname)

        return self.change_equation_name(oldname, newname)

class CellMLParser(object):
    """
    This module parses a CellML XML-file and converts it to PyCC code
    """
    @staticmethod
    def default_parameters():
        return ParameterDict(
            change_state_names=Param([], description="A list of state names "\
                                     "which should be changed to not interfere "\
                                     "with for example a parameter name."), 
            grouping=OptionParam("encapsulation", ["encapsulation", "containment"], \
                                 description="Determines what type of grouping "\
                                 "should be used when the cellml model is parsed."),
            use_sympy_integers=Param(False, description="If yes dedicated sympy "\
                                     "integers will be used instead of the "\
                                     "integers 0-10. This will turn 1/2 into a "\
                                     "sympy rational instead of the float 0.5."),
            )
    
    def __init__(self, model_source, targets=None, params=None):
        """
        Arguments:
        ----------
        
        model_source: str
            Path or url to CellML file
        targets : list, dict (optional)
            Components of the model to parse
        params : dict
            A dict with parameters for the 
        """

        targets = targets or []
        params = params or {}
        check_arg(model_source, str)
        check_arg(targets, (list, dict))
        self._params = self.default_parameters()
        self._params.update(params)

        # Open file or url
        try:
            fp = open(model_source, "r")
        except IOError:
            try:
                fp = urllib.urlopen(model_source)
            except:
                error("ERROR: Unable to open " + model_source)

        self.model_source = model_source
        self.cellml = ElementTree.parse(fp).getroot()
        self.mathmlparser = MathMLBaseParser(self._params.use_sympy_integers)
        self.cellml_namespace = self.cellml.tag.split("}")[0] + "}"
        self.name = self.cellml.attrib['name']
        self.parse_units()
        self.components = self.parse_components(targets)
        self.documentation = self.parse_documentation()

    def parse_documentation(self):
        """
        Parse the documentation of the article
        """
        namespace = "{http://cellml.org/tmp-documentation}"
        article = self.cellml.getiterator(namespace+"article")
        if not article:
            return ""

        article = article[0]

        # Get title
        if article.getiterator(namespace+"articleinfo") and \
               article.getiterator(namespace+"articleinfo")[0].\
               getiterator(namespace+"title"):
            title = article.getiterator(namespace+"articleinfo")[0].\
                    getiterator(namespace+"title")[0].text
        else:
            title = ""

        # Get model structure comments
        for child in article.getchildren():
            if child.attrib.get("id") == "sec_structure":
                content = []
                for par in child.getiterator(namespace+"para"):
                    # Get lines
                    splitted_line = deque(("".join(text.strip() \
                                        for text in par.itertext())).\
                                          split(" "))

                    # Cut them in lines which are not longer than 80 characters
                    ret_lines = []
                    while splitted_line:
                        line_stumps = []
                        line_length = 0 
                        while splitted_line and (line_length + \
                                                 len(splitted_line[0]) < 80):
                            line_stumps.append(splitted_line.popleft())
                            line_length += len(line_stumps[-1]) + 1
                        ret_lines.append(" ".join(line_stumps))
                    
                    content.extend(ret_lines)
                    content.append("\n")

                # Clean up content
                content = ("\n".join(cont.strip() for cont in content)).\
                          replace("  ", " ").replace(" .", ".").replace(\
                    " ,", ",")
                break
        else:
            content = ""
            
        if title or content:
            return "%s\n\n%s" % (title, content)

        return ""


    def _gettag(self, node):
        """
        Splits off the namespace part from name, and returns the rest, the tag
        """
        return "".join(node.tag.split("}")[1:])
    
    def parse_units(self):
        """
        Parse any declared units in the model
        """
        collected_units = OrderedDict()
        unit_map = si_unit_map.copy()

        # Extend unit_map
        for cellml_pref, pref in prefix_map.items():
            if cellml_pref:
                unit_map[cellml_pref+"molar"] = pref+"M"
                collected_units[cellml_pref+"molar"] = {pref+"M": (pref+"M", "1")}
        
        for units in self.get_iterator("units"):
            unit_name = units.attrib["name"]
            
            if unit_name in unit_map:
                continue

            collected_parts = OrderedDict()
            for unit in units.getchildren():
                if unit.attrib.get("multiplier"):
                    warning("skipping multiplier in unit {0}".format(units.attrib["name"]))
                if unit.attrib.get("multiplier"):
                    warning("skipping multiplier in unit {0}".format(units.attrib["name"]))
                cellml_unit = unit.attrib.get("units")
                prefix = prefix_map[unit.attrib.get("prefix")]
                exponent = unit.attrib.get("exponent", "1")
                if cellml_unit in si_unit_map:
                    abbrev = si_unit_map[cellml_unit]
                    name = prefix+abbrev
                    if exponent not in ["0", "1"]:
                        fullname = name + "**" + exponent
                    else:
                        fullname = name
            
                    collected_parts[name] = (fullname, exponent)
                elif cellml_unit in collected_units:
                    if prefix:
                        warning("Skipping prefix of unit '{0}'".format(cellml_unit))
                    for name, (fullnam, part_exponent) in collected_units[cellml_unit].items():
                        new_exponent = str(int(part_exponent) * int(exponent))
                        if new_exponent not in ["0", "1"]:
                            fullname = name + "**" + new_exponent
                        else:
                            fullname = name
                        
                        collected_parts[name] = (fullname, exponent)
                    
                else:
                    warning("Unknown unit '{0}'".format(cellml_unit))
                    break
            
            collected_units[unit_name] = collected_parts
            unit_map[unit_name] = "*".join(fullname for fullname, exp \
                                           in collected_parts.values())
           
        # Store unit map
        self.units_map = unit_map

    def get_iterator(self, name, item=None):
        """
        Return an element tree iterator 
        """
        
        item = item if item is not None else self.cellml
        return item.getiterator(self.cellml_namespace+name)

    def check_and_register_component_variables(\
        self, comp, collected_states, collected_parameters, collected_equations):
        """
        Check if component variables are allready collected
        """
        
        # Check for duplication of states
        for name in comp.state_variables.keys():
            der_name = "d{0}_dt".format(name)
            if name in self._params.change_state_names:
                newname = name + "_" + comp_name.split("_")[0]
                comp.change_state_name(name, newname)
                name = newname

            # Check state name vs collected state names
            elif name in collected_states:
                state_comp = collected_states[name]
                warning("Same state name: '{0}' is used in component: '{1}' "\
                        "and '{2}'.".format(name, comp.name, state_comp.name))
                for change_comp in [comp, state_comp]:
                    if change_comp.state_variables[name]["private"] and \
                           change_comp.equation_info[der_name]["private"]:
                        new_name = change_comp.change_state_name(name)
                        if change_comp == state_comp:
                            collected_states[new_name] = change_comp
                        break
                else:
                    warning("Could not resolve duplicated state name {0} in "\
                            "component {1} and {2}. None of them are private "\
                            "to the components.".format(name, comp.name,
                                                      state_comp.name))

            # Check state name vs collected parameter names
            elif name in collected_parameters:
                param_comp = collected_parameters[name]
                warning("State name: '{0}' from component '{1}' is used as "\
                        "parameter in component '{2}'.".format(\
                            name, comp.name, param_comp.name))

                # If parameter is private we change that
                if param_comp.parameters[name]["private"]:
                    new_name = param_comp.change_parameter_name(name)
                    collected_parameters.pop(name)
                    collected_parameters[new_name] = param_comp

                # Elseif the state variable is private we change that one
                elif comp.state_variables[name]["private"] and \
                         comp.equation_info[der_name]["private"]:
                    name = comp.change_state_name(name)

                else:
                    warning("Could not resolve duplicated state and "\
                            "parameter name {0} in component {1} and {2}. "\
                            "None of them are private to the "\
                            "components.".format(name, comp.name,
                                               param_comp.name))
                            
            # Check state name vs collected equation names
            elif name in collected_equations:
                eq_comp = collected_equations[name]
                warning("State name '{0}' from component '{1}' is used as "\
                        "parameter in component '{2}'.".format(\
                            name, comp.name, eq_comp.name))

                # If state is private we change it
                if comp.state_variables[name]["private"] and \
                       comp.equation_info[der_name]["private"]:
                    name = comp.change_state_name(name)
                elif eq_comp.equation_info[name]["private"]:
                    new_name = eq_comp.change_equation_name(name)
                    collected_equation.pop(name)
                    collected_equation[new_name] = eq_comp

                else:
                    warning("Could not resolve duplicated state and "\
                            "equation name {0} in component {1} and {2}. "\
                            "None of them are private to the "\
                            "components.".format(name, comp.name, eq_comp.name))
                            
            # Register the collected state
            collected_states[name] = comp

        # Check parameters
        for name in comp.parameters.keys():

            # Check parameter name vs collected state names
            if name in collected_states:
                state_comp = collected_states[name]
                warning("Same parameter and state name: '{0}' is used in "\
                        "component '{1}' and '{2}'.".format(\
                            name, comp.name, state_comp.name))
                # If parameter is private we change that
                if comp.parameters[name]["private"]:
                    name = comp.change_parameter_name(name)
                elif state_comp.state_variables[name]["private"] and \
                         state_comp.equation_info[der_name]["private"]:
                    new_name = state_comp.change_state_name(name)
                    collected_states.pop(name)
                    collected_states[new_name] = state_comp
                    
                else:
                    warning("Could not resolve duplicated state name {0} in "\
                            "component {1} and {2}. None of them are private "\
                            "to the components.".format(name, comp.name,
                                                      state_comp.name))

            # Check state name vs collected parameter names
            elif name in collected_parameters:
                param_comp = collected_parameters[name]
                warning("Parameter name '{0}' from component '{1}' is used as "\
                        "parameter in component '{2}'.".format(\
                            name, comp.name, param_comp.name))

                # If registered parameter is private we change that
                if param_comp.parameters[name]["private"]:
                    new_name = param_comp.change_parameter_name(name)
                    collected_parameters.pop(name)
                    collected_parameters[new_name] = param_comp

                # if the parameter is private we change that one
                elif comp.parameters[name]["private"]:
                    name = comp.change_parameter_name(name)

                else:
                    warning("Could not resolve duplicated parameter names "\
                            "{0} in component {1} and {2}. "\
                            "None of them are private to the "\
                            "components.".format(name, comp.name,
                                               param_comp.name))
                            
            # Check state name vs collected equation names
            elif name in collected_equations:
                eq_comp = collected_equations[name]
                warning("Parameter name '{0}' from component '{1}' "\
                        "is used as parameter in component '{2}'.".format(\
                            name, comp.name, eq_comp.name))

                # If parameter is private we change it
                if comp.parameters[name]["private"]:
                    name = comp.change_parameter_name(name)
                elif eq_comp.equation_info[name]["private"]:
                    new_name = eq_comp.change_equation_name(name)
                    collected_equation.pop(name)
                    collected_equation[new_name] = eq_comp

                else:
                    warning("Could not resolve duplicated parameter and "\
                            "equation names {0} in component {1} and {2}. "\
                            "None of them are private to the "\
                            "components.".format(name, comp.name, eq_comp.name))
                            
            collected_parameters[name] = comp

        # Check equation name
        for eq in comp.equations:
            name = eq.name

            # If the equation is private we do not care if it is duplicated
            # elsewhere. Risky...
            if comp.equation_info[name]["private"]:
                collected_equations[name] = comp
                continue
            
            # Check equation name vs collected state names
            if name in collected_states:
                state_comp = collected_states[name]
                warning("Same equation and state name '{0}' is used in "\
                        "component '{1}' and '{2}'.".format(\
                            name, comp.name, state_comp.name))
                # If equation is private we change that
                if comp.equation_info[name]["private"]:
                    name = comp.change_equation_name(name)
                elif state_comp.state_variables[name]["private"] and \
                         state_comp.equation_info[der_name]["private"]:
                    new_name = state_comp.change_state_name(name)
                    collected_states.pop(name)
                    collected_states[new_name] = state_comp
                    
                else:
                    warning("Could not resolve duplicated state and "\
                            "equation name {0} in component {1} and {2}. "\
                            "None of them are private to the "\
                            "components.".format(name, comp.name,
                                               state_comp.name))

            # Check state name vs collected parameter names
            elif name in collected_parameters:
                param_comp = collected_parameters[name]
                warning("Equation name '{0}' from component '{1}' is used as "\
                        "parameter in component '{2}'." % (\
                            name, comp.name, param_comp.name))

                # If registered parameter is private we change that
                if param_comp.parameters[name]["private"]:
                    new_name = param_comp.change_parameter_name(name)
                    collected_parameters.pop(name)
                    collected_parameters[new_name] = param_comp

                # If the parameter is private we change that one
                elif comp.equation_info[name]["private"]:
                    name = comp.change_equation_name(name)

                else:
                    warning("Could not resolve duplicated parameter and "\
                            "equation name {0} in component {1} and {2}. "\
                            "None of them are private to the "\
                            "components.".format(name, comp.name,
                                               param_comp.name))
                            
            # Check equation name vs collected equation names
            elif name in collected_equations:
                eq_comp = collected_equations[name]
                warning("Equation name '{0}' from component '{1}' is used as "\
                        "parameter in component '{2}'.".format(\
                            name, comp.name, eq_comp.name))

                # If parameter is private we change it
                if comp.equation_info[name]["private"]:
                    name = comp.change_equation_name(name)
                    
                elif eq_comp.equation_info[name]["private"]:
                    new_name = eq_comp.change_equation_name(name)
                    collected_equation.pop(name)
                    collected_equation[new_name] = eq_comp

                else:
                    warning("Could not resolve duplicated parameter and "\
                            "equation names {0} in component {1} and {2}. "\
                            "None of them are private to the "\
                            "components.".format(name, comp.name, eq_comp.name))
            
            collected_equations[name] = comp
            
    def parse_imported_model(self):
        """
        Parse any imported models
        """

        components = OrderedDict()

        # Collect states and parameters to check for duplicates
        collected_states = dict()
        collected_parameters = dict()
        collected_equations = dict()

        # Import other models
        for model in self.get_iterator("import"):
            import_comp_names = dict()

            for comp in self.get_iterator("component", model):
                
                import_comp_names[comp.attrib["component_ref"]] = \
                                                        comp.attrib["name"]

            model_parser = CellMLParser(\
                model.attrib["{http://www.w3.org/1999/xlink}href"], \
                import_comp_names)

            for comp in model_parser.components:
                components[comp.name] = comp

        # Extract states, parameters and equations
        for comp in components.values():

            self.check_and_register_component_variables(\
                comp, collected_states, collected_parameters, \
                collected_equations)
            
        return components, collected_states, collected_parameters,\
               collected_equations

    def get_parents(self, grouping, element=None):
        """
        If group was used in the cellml use it to gather parent information
        about the components
        """
        # Collect component encapsulation
        def get_encapsulation(elements, all_parents, parent=None):
            children = {}
            for encap in elements:
                name = encap.attrib["component"]
                all_parents[name] = parent
                if encap.getchildren():
                    nested_children = get_encapsulation(\
                        encap.getchildren(), all_parents, name)
                    children[name] = dict(children=nested_children, \
                                          parent=parent)
                else:
                    children[name] = dict(children=None, parent=parent)

            return children

        encapsulations = dict()
        all_parents = dict()
        for group in self.get_iterator("group", element):
            children = group.getchildren()

            if children and children[0].attrib.get("relationship") == \
                   grouping:
                encapsulations = get_encapsulation(children[1:], all_parents)

        # If no group information in cellml extract potential parent information
        # from component names
        if not all_parents:

            # Iterate over the components
            comp_names = [comp.attrib["name"] for comp in self.get_iterator(\
                "component", element)]

            for parent_name in comp_names:
                for name in comp_names:
                    if parent_name in name and parent_name != name:
                        all_parents[name] = parent_name

        return encapsulations, all_parents

    def parse_single_component(self, comp, collected_states, \
                               collected_parameters, collected_equations):
        """
        Parse a single component and create a Component object
        """
        comp_name = comp.attrib["name"]
        
        # Collect variables and equations
        variables = OrderedDict()
        equations = []
        state_variables = OrderedDict()
        #derivatives = []

        # Get variable and initial values
        for var in self.get_iterator("variable", comp):

            var_name = var.attrib["name"]
            if var_name in _all_keywords:
                var_name = var_name + "_"

            # Store variables using initial and unit
            variables[var_name] = dict(\
                init=var.attrib.get("initial_value"),\
                unit=self.units_map[var.attrib.get("units")],\
                private=not (var.attrib.get("public_interface")=="out" or \
                             var.attrib.get("private_interface")=="out"))

        # Get equations
        for math in comp.getiterator(\
            "{http://www.w3.org/1998/Math/MathML}math"):
            for eq in math.getchildren():
                equation_list, state_variable, derivative, \
                        used_variables = self.mathmlparser.parse(eq)

                # Get equation name
                eq_name = equation_list[0]
                
                if eq_name in _all_keywords:
                    equation_list[0] = eq_name + "_"
                    eq_name = equation_list[0]
                
                # Discard collected equation name from used variables
                used_variables.discard(eq_name)
                
                assert(re.findall("(\w+)", eq_name)[0]==eq_name)
                assert(equation_list[1] == self.mathmlparser["eq"])
                equations.append(Equation(eq_name, equation_list[2:],
                                          used_variables))
                
                # Do not register state variables twice
                if state_variable is not None and \
                       state_variable not in state_variables:
                    state_variables[state_variable] = derivative

        # Create Component
        comp = Component(comp_name, variables, equations, state_variables)

        # Collect and check variables
        self.check_and_register_component_variables(\
            comp, collected_states, collected_parameters, \
            collected_equations)

        return comp

    def sort_components(self, components):
        
        # Check internal dependencies
        for comp0 in components:
            for comp1 in components:
                if comp0 == comp1:
                    continue
                comp0.check_dependencies(comp1)

        def simple_sort(components):
            components = deque(components)
            dependant_components = []
            sorted_components = []
            while components:
                component = components.popleft()

                # Chek for circular dependancy
                if dependant_components.count(component) > 4:
                    components.append(component)
                    break
                
                if any(dep in components for dep in \
                       component.dependencies):
                    components.append(component)
                    dependant_components.append(component)
                else:
                    sorted_components.append(component)

            return sorted_components, list(components)

        # Initial sorting
        sorted_components, circular_components = simple_sort(components)

        # If no circular dependencies
        if not circular_components:
            return sorted_components

        try:
            import networkx as nx
        except:
            warning("networkx could not be imported. Circular "\
                    "dependencies between components will not be sorted.")
            return sorted_components + circular_components
        
        # Collect zero and one dependencies
        zero_dep_equations = set()
        one_dep_equations = set()
        equation_map = {}
        
        # Gather zero dependent equations
        for comp in circular_components:
            for dep_comp, equations in comp.dependencies.items():
                #if dep_comp not in circular_components:
                #    continue
                for equation in equations:
                    if not equation.dependent_equations and equation.name \
                           in comp.used_variables:
                        zero_dep_equations.add(equation)
                        equation_map[equation.name] = equation

        # Check for one dependency if that is the zero one
        one_dep_zero_dep = defaultdict(set)
        for comp in circular_components:
            for dep_comp, equations in comp.dependencies.items():
                #if dep_comp not in circular_components:
                #    continue
                for equation in equations:
                    if len(equation.dependent_equations) == 1 and \
                           equation.name in comp.used_variables and \
                           equation.dependent_equations[0] in \
                           zero_dep_equations:
                        equation_map[equation.name] = equation
                        one_dep_equations.add(equation)
                        one_dep_zero_dep[equation.name].add(\
                            equation.dependent_equations[0].name)

        # Try to eliminate circular dependency
        # Extract dependent equation to a new component
        ode_comp = Component(self.name, {}, [], {})

        # Valid edges for removal
        valid_edges = [eq.name for eq in zero_dep_equations] + \
                      [eq.name for eq in one_dep_equations]
        
        G = nx.MultiDiGraph()
        G.add_nodes_from([comp.name for comp in components])
        
        # Build graph
        for comp in components:
            [G.add_edge(dep.name, comp.name, key=equation.name)
             for dep, equations in comp.dependencies.items() \
             for equation in equations]

        # collecte edges that breaks cycles
        cycle_breakers = []

        # Collect data over the best edge to remove
        edge_score = defaultdict(lambda : 0)
        edge_to_nodes = defaultdict(set)
        for cycle in nx.simple_cycles(G):
            local_breaker = []
            for n0, n1 in zip(cycle[:-1], cycle[1:]):
                if all(edge in valid_edges for edge in G[n0][n1]):
                    local_breaker.append([edge for edge in G[n0][n1]])
                    for edge in local_breaker[-1]:
                        edge_to_nodes[edge].add((n0, n1))
                        
            cycle_breakers.append(local_breaker)
            for local_edges in local_breaker:
                for edge in local_edges:
                    edge_score[edge] += 1

        # Sort the edges we should remove by a score given by how many
        # cycles it breaks by being removed.
        edge_score = sorted((edge_score[edge], edge) for edge in edge_score)

        # Collect edges to be removed
        edge_removal = set()
        cycles_fixed = [False]*len(cycle_breakers)
        while not all(cycles_fixed):
            score, edge = edge_score.pop()

            # Remove this/these edge/edges this iteration
            local_removal = [edge]

            # If adding a one dep edge we need to also remove its dependent edge
            if edge in one_dep_zero_dep:
                local_removal.extend(one_dep_zero_dep[edge])

            for edge_remove in local_removal:
                for i, local_breakers in enumerate(cycle_breakers):

                    # If the cycle is fixed
                    if cycles_fixed[i]:
                        continue

                    # Go through the collected valid breakers
                    for j, local_edges in enumerate(local_breakers):

                        # If the removed edge in in the local edges
                        if edge_remove in local_edges:
                            local_edges.remove(edge_remove)
                            if len(local_edges) == 0:
                                cycles_fixed[i] = True
                                edge_removal.add(edge_remove)
                                break

        # Remove the edges from the graph
        for edge in edge_removal:
            for n0, n1 in edge_to_nodes[edge]:
                G.remove_edge(n0, n1, key=edge)

        # Move the marked edges (equations) from the relevant components
        removed_equations = {}
        for edge in edge_removal:

            eq = equation_map[edge]
            
            old_comp = eq.component
            ode_comp.equations.append(eq)
            old_comp.equations.remove(eq)

            # Store changed 
            removed_equations[eq] = old_comp

            # Transfer used_in to new component
            new_dependent_componets = OrderedDict()
            for dep_comp, equations in old_comp.used_in.items():
                assert equations

                if eq not in equations or dep_comp == ode_comp:
                    continue

                # Remove dependency from old component and add it to the new
                if len(equations) == 1:
                    new_dependent_componets[dep_comp] = \
                                    old_comp.used_in.pop(dep_comp)
                else:
                    new_dependent_componets[dep_comp] = [\
                        old_comp.used_in[dep_comp].pop(equations.index(eq))]

                # Change component dependencies
                if old_comp in dep_comp.dependencies and eq in \
                       dep_comp.dependencies[old_comp]:
                    if len(dep_comp.dependencies[old_comp]) == 1:
                        dep_comp.dependencies.pop(old_comp)
                    else:
                        dep_comp.dependencies[old_comp].remove(eq)

                if ode_comp not in dep_comp.dependencies:
                    dep_comp.dependencies[ode_comp] = [eq]
                else:
                    dep_comp.dependencies[ode_comp].append(eq)

            # Store new component to the extracted equation
            eq.component = ode_comp
            ode_comp.dependent_componets = new_dependent_componets
            
            if ode_comp not in sorted_components:
                sorted_components.insert(0, ode_comp)
            
        # Sort newly added equations
        ode_comp.sort_and_store_equations(ode_comp.equations)

        # Sort graph components and apply the sortings to the collected
        # CellML compoents
        sorted_components = nx.topological_sort(G)
        components.sort(lambda n0, n1: cmp(sorted_components.index(n0.name), \
                                           sorted_components.index(n1.name)))

        # Insert the ODE component with extracted equations first
        components.insert(0, ode_comp)

        warning("To avoid circular dependency the following equations "\
                "has been moved:")
        
        for eq, old_comp in removed_equations.items():
            warning("{0} : from {1} to {2} component".format(\
                eq.name, old_comp.name, ode_comp.name))

        return components

    def parse_components(self, targets):
        """
        Build a dictionary containing dictionarys describing each
        component of the cellml model
        """

        # Parse imported components
        components, collected_states, collected_parameters, \
                    collected_equations = self.parse_imported_model()

        # Get parent relationship between components
        encapsulations, all_parents = self.get_parents(self._params.grouping)

        if targets:

            # If the parent information was not of type encapsulation
            # regather parent information
            if self._params.grouping != "encapsulation":

                encapsulations, dummy = self.get_parents("encapsulation")

            # Add any encapsulated components to the target list
            for target, new_target_name in targets.items():
                if target in encapsulations:
                    for child in encapsulations[target]["children"]:
                        targets[child] = child.replace(\
                            target, new_target_name)

            target_parents = dict()

            # Update all_parents
            for comp_name, parent_name in all_parents.items():
                if parent_name not in targets:
                    continue
                target_parents[targets[comp_name]] = targets[parent_name]
        
        # Iterate over the components
        for comp in self.get_iterator("component"):
            comp_name = comp.attrib["name"]

            # Only parse selected and non-empty components
            if (targets and comp_name not in targets) or \
                   len(comp.getchildren()) == 0:
                continue

            # If targets provides a name mapping give the component a new name
            if targets and isinstance(targets, dict):
                new_name = targets[comp_name]
                comp.attrib["name"] = new_name
                comp_name = new_name
            
            # Store component
            components[comp_name] = self.parse_single_component(\
                comp, collected_states, collected_parameters, collected_equations)

        # Add parent information
        for name, comp in components.items():

            if targets:
                parent_name = target_parents.get(name)
            else:
                parent_name = all_parents.get(name)

            if parent_name:
                comp.parent = components[parent_name]
                
                # If parent name in child name, reduce child name length
                if parent_name in comp.name:
                    comp.name = comp.name.replace(parent_name, "").strip("_")

                components[parent_name].children.append(comp)

        # If we only extract a sub set of component we do not sort
        if targets:
            return components.values()

        # Before dependencies are checked we change names according to
        # variable mappings in the original CellML file
        new_variable_names, same_variable_names = self.parse_name_mappings()
        for comp, variables in new_variable_names.items():

            # Iterate over old and new names
            for oldname, newnames in variables.items():
                
                # Check if the oldname is used in any components
                oldname_used = oldname in same_variable_names[comp]
                
                # If there are only one newname we change the name of the
                # original equation
                if len(newnames) == 1:

                    if not oldname_used:
                        newname = newnames.keys()[0]
                        if components[comp].get_variable_type(oldname) is not None:
                            components[comp].change_variable_name(oldname, newname)
                        else:
                            for child in components[comp].children:
                                if child.get_variable_type(oldname) is not None:
                                    child.change_variable_name(oldname, newname)
                                    break
                    else:
                        # FIXME: Add equation with name change to component
                        pass
                else:
                    # FIXME: Add equation with name change to component
                    pass
        
        # Add dependencies and sort the components accordingly
        return self.sort_components(components.values())
        
    def parse_name_mappings(self):
        new_variable_names = dict()
        same_variable_names = dict()
        
        for con in self.get_iterator("connection"):
            con_map = self.get_iterator("map_components", con)[0]
            comp1 = con_map.attrib["component_1"]
            comp2 = con_map.attrib["component_2"]
            
            for var_map in self.get_iterator("map_variables", con):
                var1 = var_map.attrib["variable_1"]
                var2 = var_map.attrib["variable_2"]
                
                if var1 != var2:
                    if comp1 not in new_variable_names:
                        new_variable_names[comp1] = {var1:defaultdict(list)}
                    elif var1 not in new_variable_names[comp1]:
                        new_variable_names[comp1][var1] = defaultdict(list)
                    new_variable_names[comp1][var1][var2].append(comp2)

                else:
                    if comp1 not in same_variable_names:
                        same_variable_names[comp1] = {var1:[comp2]}
                    elif var1 not in same_variable_names[comp1]:
                        same_variable_names[comp1][var1] = [comp2]
                    else:
                        same_variable_names[comp1][var1].append(comp2)

        return new_variable_names, same_variable_names
    
    def to_gotran(self):
        """
        Generate a gotran file
        """
        gotran_lines = []
        for docline in self.documentation.split("\n"):
            gotran_lines.append("# " + docline)

        if gotran_lines:
            gotran_lines.extend([""])

        # Add component info
        declaration_lines = []
        equation_lines = []

        def unders_score_replace(comp):

            new_name = comp.name.replace("_", " ")

            # If only 1 state it might be included in the name
            state_name = ""
            if len(comp.state_variables) == 1:
                state_name = comp.state_variables.keys()[0]
                if state_name in comp.name:
                    new_name = state_name.join(part.replace("_", " ") \
                                for part in comp.name.split(state_name))

            # Captitalize first word
            single_words = new_name.split(" ")
            if len(single_words[0]) > 1 and "_" not in single_words[0]:
                single_words[0] = single_words[0][0].upper()+single_words[0][1:]

            # If first word is only a character we assume the first two
            # words to stick together
            elif len(single_words[0]) == 1 and len(single_words)>1 and \
                     single_words[0] != state_name:
                single_words = [single_words[0]+"_"+single_words[1]] + single_words[2:]

            return " ".join(single_words)

        # Iterate over components and collect stuff
        for comp in self.components:
            
            names = deque([unders_score_replace(comp)])

            parent = comp.parent
            while parent is not None:
                names.appendleft(unders_score_replace(parent))
                parent = parent.parent

            comp_name = ", ".join("\"{0}\"".format(name) for name in names)
            
            # Collect initial state values
            if comp.state_variables:
                declaration_lines.append("")
                declaration_lines.append("states({0},".format(comp_name))
                for name, info in comp.state_variables.items():
                    if info["unit"] != "1":
                        declaration_lines.append("       {0} = ScalarParam({1}"\
                                ", unit=\"{2}\"),".format(name, info["init"], info["unit"]))
                    else:
                        declaration_lines.append("       {0} = {1},".format(\
                            name, info["init"]))
                declaration_lines[-1] = declaration_lines[-1][:-1]+")"

            # Collect initial parameters values
            if comp.parameters:
                declaration_lines.append("")
                declaration_lines.append("parameters({0},".format(comp_name))
                for name, info in comp.parameters.items():
                    if info["unit"] != "1":
                        declaration_lines.append("           {0} = ScalarParam({1}"\
                                ", unit=\"{2}\"),".format(name, info["init"], info["unit"]))
                    else:
                        declaration_lines.append("           {0} = {1},".format(\
                            name, info["init"]))
                declaration_lines[-1] = declaration_lines[-1][:-1]+")"

            # Collect all intermediate equations
            if comp.equations:
                equation_lines.append("")
                equation_lines.append("component({0})".format(\
                    comp_name))

                equation_lines.extend("{0} = {1}".format(eq.name, "".join(eq.expr))\
                                      for eq in comp.equations)

        gotran_lines.append("# gotran file generated by cellml2gotran from "\
                            "{0}".format(self.model_source))
        gotran_lines.extend(declaration_lines)
        gotran_lines.extend(equation_lines)
        gotran_lines.append("")
        gotran_lines.append("")
        

        # Return joined lines
        return "\n".join(gotran_lines)

        # Write file
        open("{0}.ode".format(self.name), \
             "w").write()

def cellml2ode(model_source, **options):
    """
    Convert a CellML model into an ode

    Arguments:
    ----------
    model_source: str
        Path or url to CellML file
    options : dict
        Optional parameters to control cellml parser
    """
    check_arg(model_source, str)
    from gotran import exec_ode
    params = CellMLParser.default_parameters()
    params.update(params)
    cellml = CellMLParser(model_source, params=params)
    return exec_ode(cellml.to_gotran(), cellml.name)
    
