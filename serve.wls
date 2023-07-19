Import["https://raw.githubusercontent.com/JerryI/tinyweb-mathematica/master/Tinyweb/Tinyweb.wl"]

server = WEBServer["addr" -> "127.0.0.1:5699", "path" -> {FileNameJoin[{"tests", "public"}], FileNameJoin[{"tests", "Packages"}], "dist"}, "socket-close"->True]
server = server // WEBServerStart
StringTemplate["Open your browser at http://``"][server["addr"]]//Print

While[True, Pause[0.1]];