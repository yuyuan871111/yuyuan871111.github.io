# Personal Academic Website

[![Deploy Hugo site to Pages](https://github.com/yuyuan871111/yuyuan871111.github.io/actions/workflows/hugo.yml/badge.svg)](https://github.com/yuyuan871111/yuyuan871111.github.io/actions/workflows/hugo.yml)

Academic portfolio and research showcase website built with [Hugo](https://gohugo.io/) and the [hugo-profile](https://github.com/gurusabarish/hugo-profile) theme, extended with custom interactive directory layouts for publications, conferences, and honors.

---

## Quick Start & Local Development

```bash
# 1. Preview locally with live reloading
hugo server -D

# 2. Open in browser
open http://localhost:1313/

# 3. Deploy changes to GitHub Pages
git add .
git commit -m "Update publications and honors"
git push
```

---

## Content Architecture Overview

| Record Type | Primary Data Location | Target Page / URL |
| :--- | :--- | :--- |
| **Publications** | `hugo.yaml` (`params.achievements.publications`) | `/blogs/publication` & Homepage projects |
| **Presentations (Oral / Poster)** | `hugo.yaml` (`oralPresentations` / `posterPresentations`) | `/blogs/publication` & Homepage projects |
| **Conferences (Attendance)** | `hugo.yaml` (`attendanceOnly`) | `/blogs/publication` |
| **Service & Teaching** | `hugo.yaml` (`service`) | `/blogs/publication` & Homepage collaborative section |
| **Honors & Awards** | `hugo.yaml` (`params.achievements.honors`) | `/blogs/honors` |
| **Projects (Homepage Showcase)** | `hugo.yaml` (`params.projects.items`) | Homepage (`/#projects`) |
| **News & Blog Posts** | `content/blogs/*.md` | Homepage Recent News & `/blogs` |
| **Images & Figures** | `static/images/blog_pics/` | Used across cards & articles |

---

## How to Update Your Records

### 1. Adding a New Publication

When you publish a new paper or preprint:

1. Open `hugo.yaml` and add an entry under `params.achievements.publications`:

```yaml
      - title: "Your Paper Title Here"
        authors: "Yang, Y. Y., Coauthor, A., and Advisor, B."
        venue: "Journal Name (or bioRxiv)"
        date: "2026"                # Or month: "Mar 2026"
        type: "Research"            # e.g., Research, Preprint, Review, Tool
        highlight: true             # Set true to feature on homepage project card
        projectTag: "DLPROBSEG"     # Link to project tag (optional)
        projectUrl: "/tags/dlprobseg"
        doi: "10.1021/acs.jcim.xxxx"
        doiUrl: "https://doi.org/10.1021/acs.jcim.xxxx"
```

2. *(Optional)* Add a formatted citation row to `content/blogs/publication.md` under `## Publication` as a markdown fallback:
```markdown
8\. [<o>research</o>, [DLPROBSEG](/tags/dlprobseg)] **Yang, Y. Y.**, ... (2026). "Your Paper Title." *Journal*, DOI: [10.1021/xxx](https://doi.org/xxx)
```

> [!TIP]
> Setting `projectTag` automatically embeds this publication inside the corresponding project card on the homepage under **Related publications**.

---

### 2. Adding a Conference Presentation or Attendance

#### A. Oral Presentations
In `hugo.yaml` under `params.achievements.oralPresentations`:
```yaml
      - title: "Talk Title Here"
        authors: "Yang, Y. Y., Advisor, A."
        event: "Conference Name 2026"
        location: "London, UK"      # Or "Virtual"
        date: "Jun 2026"
        type: "Oral"
        highlight: true             # Shows in project card if projectTag matches
        projectTag: "DLPROBSEG"     # Optional
        projectUrl: "/tags/dlprobseg"
        eventUrl: "https://..."     # Optional link to event / recording
```

#### B. Poster Presentations
In `hugo.yaml` under `params.achievements.posterPresentations`:
```yaml
      - title: "Poster Title Here"
        authors: "Yang, Y. Y., Advisor, A."
        event: "Symposium Name 2026"
        location: "Oxford, UK"
        date: "Nov 2026"
        type: "Poster"              # Or "Poster / Lightning Talk"
        highlight: true
        projectTag: "DLPROBSEG"
        eventUrl: "https://..."
```

#### C. Conference Attendance Only
In `hugo.yaml` under `params.achievements.attendanceOnly`:
```yaml
      - title: "Conference Name 2026"
        event: "Organizing Body or Institution"
        location: "London, UK"
        date: "Oct 2026"
        year: "2026"
        type: "Attendance"
        url: "https://..."          # Optional
```

---

### 3. Adding an Honor, Grant, or Award

Open `hugo.yaml` and add an entry under `params.achievements.honors`:

```yaml
      - title: "Award or Scholarship Name"
        organization: "Granting Institution / Body"
        program: "Specific Programme or Department"
        date: "2026"
        year: "2026"
        category: "grant"           # Options: grant, award, competition, distinction
        categoryLabel: "Grants & Fellowships"  # Matches filter pill
        badge: "PhD Studentship"    # Distinction pill text (e.g. 1st Place, Best Thesis)
        location: "London, UK"
        icon: "fas fa-trophy"       # e.g., fas fa-trophy, fas fa-award, fas fa-medal, fas fa-graduation-cap
        content: "One or two sentences describing the honor, selectivity, or funded scope."
        url: "https://..."          # Optional link to certificate / DOI / announcement
        urlText: "View Announcement" # Optional button text
```

#### Supported Honor Categories:
* `grant` → **Grants & Fellowships** (icon: `fas fa-graduation-cap` or `fas fa-flask`)
* `award` → **Awards & Thesis** (icon: `fas fa-trophy`)
* `competition` → **Competitions & Hackathons** (icon: `fas fa-award` or `fas fa-medal`)
* `distinction` → **Academic Distinctions** (icon: `fas fa-star`)

---

### 4. Adding a Research Project to the Homepage

In `hugo.yaml` under `params.projects.items`:

```yaml
      - title: "Project Name Here"
        projectTag: "TAGNAME"        # Used for linking publications and talks
        badges:
          - "Methodology"
          - "Tool / Package"
        content: |-
          One concise paragraph explaining the scientific motivation, methodology, and key finding.
        image: /images/blog_pics/your_study_figure.png
        links: 
          - icon: fa-solid fa-blog
            url: /tags/tagname
          - icon: fab fa-github
            url: https://github.com/yuyuan871111/repo_name
        featured:                    # Optional live web service button
          name: "Live Web App"
          link: "https://..."
```

> [!NOTE]
> Project figures are displayed in an uncropped landscape container (`180px` height). For best results, place images in `static/images/blog_pics/` and reference them as `/images/blog_pics/<filename>`.

---

### 5. Writing a News / Blog Post

Create a new markdown file in `content/blogs/<post-slug>.md`:

```markdown
---
title: "Article or News Headline"
date: 2026-03-01T10:00:00Z
draft: false
author: "Yu-Yuan (Stuart) Yang"
tags:
  - TAGNAME
image: /images/blog_pics/your_image.png
description: "Brief summary shown on listings and preview cards."
toc: true
---

Your content in Markdown goes here.
```

* Regular posts in `content/blogs/` automatically appear in the **Recent News** feed on the homepage and in the `/blogs` index.
* `publication.md` and `honors.md` have `build: { list: never }` in their frontmatter, keeping them accessible at their standalone URLs (`/blogs/publication` and `/blogs/honors`) while hiding them completely from the `/blogs` article list and pagination.

---

## Directory Reference

```text
├── content/
│   ├── blogs/
│   │   ├── publication.md     # Publications & conferences page entry
│   │   ├── honors.md          # Honors & awards page entry
│   │   └── *.md               # News articles & project writeups
├── layouts/
│   ├── _default/
│   │   ├── publications.html  # Search & filter layout for publications
│   │   └── honors.html        # Search & filter layout for honors & awards
│   └── partials/
│       ├── head.html          # Global theme styling overrides
│       └── sections/
│           ├── projects.html  # Homepage project cards + collaborative block
│           ├── education.html # Education section
│           └── footer/
│               └── recentBlogPosts.html # Homepage Recent News feed
├── static/
│   ├── images/
│   │   └── blog_pics/         # Figures, study diagrams, and card graphics
│   └── files/                 # CV and downloadable PDF documents
└── hugo.yaml                  # Central site configuration and record database
```
