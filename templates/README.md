# for developers

We are using the [`format` method](https://docs.python.org/3/tutorial/inputoutput.html#the-string-format-method) to template these files
.
So in latex you will need to write `\frac{{num}}{{den}}`  instead of `\frac{num}{den}` to make python ignore tose keywords

This should help you to visualize a better form for the Markdown and import it to be compiled with a simple procedure:

```python
with open(filename, "r") as f:
    inp = f.read()
output={"key1":"text"}
output=["key2"]:"text"
"""things happen here"""
# output like:
with open("output.md", "w+") as f:
    f.write(inp.format(**output))
```