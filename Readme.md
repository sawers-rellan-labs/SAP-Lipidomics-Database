Use this project template to organize your projects. Inspired by work of many others. Check:

- [NiceR Code](https://nicercode.github.io/blog/2013-04-05-projects/)
- [Carl Boettiger](http://www.carlboettiger.info/2012/05/06/research-workflow.html)
- [William Stafford Noble´s guide on project organization](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000424)

To start using this template download the repository from [here](https://www.dropbox.com/sh/8aurm6lbg8xc76r/AABdrJEGX0Ls08EFG0nivfRga?dl=0), change repository name and start doing great science

Some tips and recommendations.

### `data`

- In data you should store only raw data.
- This data will be read when you do your analysis but can never be modified.
- If you want to create modified data-sets from this data save it in the output folder.
- Create a document where you explain how you obtained the data and also the meaning of all your variables.

### `docs`

Here is where your article/thesis goes. You can use Latex, or Markdown. If your article is heavy on R analysis you could use Rmarkdown to directly write your article then [`knitcitations`](https://github.com/cboettig/knitcitations) is highly recommended to handle citations.  
My personal workflow involves using Atom locally, write directly in Github´s editor or using [Prose.io](prose.io). I use magic citations (that work on both Atom and online with Chrome) from [Papers](papersapp.com) to insert citations. I will be happy to buy you a student Papers license if you want to go that route. To generate References and export to Latex/Pdf via [`pandoc`](http://pandoc.org/) I use a small `makefile` script. Please check Karl Broman´s tutorial on [`make`](http://kbroman.org/minimal_make/). I have included a sample in the docs folder. With this strategy you:
- Write in Markdown (you can work locally or online) which I think is easier for beginners.
- Insert citations with Magic Citations from Papers. But you can use Mendeley or other solutions. Is easy to generate a bibtex bibliography file in Papers or Mendeley and then copy it in a bibliography.bib file.  
- Keep your article version controlled by committing frequently to Github.
- Generate PDF or latex file using a simple makefile.

### `figures`

Save here your final figures. The ones that you will use to generate article PDF.

### `output`

- Output is for data you generate with your scripts, figures etc. The idea is that everything that you put here is disposable because it can be regenerated via scripts.
- You could subdivide this folder in two: figures and data.

### `scripts`

Scripts of your analysis. If in R, Rmarkdown is highly recommended. If the project is simple you could have a single Rmarkdown file stored in the root of the repository and skip this part. An `Rproject` file is also generated in the root of the repository

### Internal communication

Finally, for each project we will generate an associated [`gitter`](gitter.im) channel as a way to keep track of ideas, suggestions etc.
