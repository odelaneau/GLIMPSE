---
title: Glossary
nav_order: 11
---

# Glossary of Terms

Below are some important concepts, software, standards, and other things you might encounter in these Docs and working with SHAPEIT5. 

{% assign terms = site.glossary | sort: "title" %}
{% for t in terms %}- [{{ t.title }}](#{{ t.title | slugify }})
{% endfor %}

{% for t in terms %}
--------

## {{ t.title }}

{% if t.link %}<{{ t.link }}>{% endif %}

{{ t.content }}

{% endfor %}

