import earthkit.regrid

print(earthkit.geo.config.get("url-download-timeout"))

with earthkit.geo.config.temporary():
    earthkit.geo.config.set("url-download-timeout", 5)
    print(earthkit.geo.config.get("url-download-timeout"))

# Temporary config can also be created with arguments:
with earthkit.geo.config.temporary("url-download-timeout", 11):
    print(earthkit.geo.config.get("url-download-timeout"))
